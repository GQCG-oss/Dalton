module xcfun_host
  use precision
#ifdef VAR_XCFUN
  use xcfun
#endif
  implicit none
  type xc_functional
    character(len=1024)         :: DFTfuncString
    integer                     :: XCFUNfunctional
    integer                     :: len
    logical                     :: gga
    logical                     :: meta
    real(realk)                 :: hfweight
    type(xc_functional),pointer :: next
  end type xc_functional

  integer,save :: XCFunFunctional
  logical :: USEXCFUN,XCFUN_DOGGA,XCFUN_DOMETA
  integer,parameter :: xcfun_type_lda=1
  integer,parameter :: xcfun_type_gga=2
  integer,parameter :: xcfun_type_metagga=3
  public :: xcfun_host_init, xcfun_host_free, USEXCFUN, &
       & xcfun_gga_xc_single_eval, xcfun2_gga_xc_single_eval, &
       & xcfun_lda_xc_eval, xcfun_gga_unres_xc_single_eval, &
       & xcfun_lda_unres_xc_single_eval, xcfun3_gga_xc_single_eval,&
       & xcfun_meta_xc_single_eval, xcfun2_lda_xc_single_eval,&
       & xcfun3_lda_xc_single_eval,&
       & xcfun_host_set_order,xcfun_type_lda,xcfun_type_gga,&
       & xcfun_type_metagga, xcfun_gga_components_xc_single_eval,&
       & xcfun_host_type, xcfundftreport, xcfun_host_augment_func
  private
  integer,save                     :: num_xc_functionals = 0
  type(xc_functional),pointer,save :: xcfun_host_first
  type(xc_functional),pointer,save :: xcfun_host_last
  type(xc_functional),pointer,save :: xcfun_host_current
  contains
  subroutine xcfun_host_init(DFTfuncString,hfweight,lupri)
    implicit none
    integer                 :: lupri
    character(*),intent(in) :: DFTfuncString
    real(realk),intent(inout) :: hfweight
    !
#ifdef VAR_XCFUN
    character(len=80),pointer     :: DFTfuncStringSingle(:)
    integer                       :: Ipos,nStrings,ierr,I,ierrLDA,ierrGGA,ierrMETA,ierrHF
    !logical                       :: GGAkeyString,ErfString           deprecated
    real(realk),pointer           :: WeightSingle(:)
    logical                       :: new_func

    !Initialization
    !GGAkeyString = .FALSE.
    !ErfString = .false.

    XCFUNfunctional = xcfun_host_get_func(DFTfuncString,hfweight)

    !Generate a new functional if it does not already exist
    IF (XCFUNfunctional.EQ.-1) THEN
       CALL xcfun_host_add_functional(DFTfuncString,hfweight,lupri)
    ENDIF
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif

  end subroutine xcfun_host_init

#ifdef VAR_XCFUN
  FUNCTION xcfun_host_get_func(DFTfuncString,hfweight)
  implicit none
  character(*) :: DFTfuncString
  real(realk)  :: hfweight
  integer      :: xcfun_host_get_func
  !
  integer :: ifunc,lenght
  type(xc_functional),pointer   :: current

  lenght = len(DFTfuncString)
  xcfun_host_get_func = -1
  current => xcfun_host_first
  DO ifunc=1,num_xc_functionals
    IF (DFTfuncString(1:lenght) == current%DFTfuncString(1:lenght)) THEN
      xcfun_host_get_func = current%XCFUNfunctional
      hfweight            = current%hfweight
      XCFUN_DOGGA         = current%gga
      XCFUN_DOMETA        = current%meta
      xcfun_host_current => current
      exit
    ENDIF
    current => current%next
  ENDDO
  END FUNCTION xcfun_host_get_func

  subroutine xcfun_host_add_functional(DFTfuncString,hfweight,lupri)
    implicit none
    integer                       :: lupri
    character(*),intent(in)       :: DFTfuncString
    real(realk),intent(inout)     :: hfweight
    !
    character(len=80),pointer     :: DFTfuncStringSingle(:)
    integer                       :: Ipos,nStrings,ierr,I,ierrLDA,ierrGGA,ierrMETA,ierrHF
    !logical                       :: GGAkeyString,ErfString           deprecated
    real(realk),pointer           :: WeightSingle(:)
    logical                       :: new_func,admm_ggacorr_single,admm_ggacorr
    type(xc_functional),pointer   :: current

    nullify(current)
    allocate(current)
    IF (num_xc_functionals.EQ.0) THEN
      xcfun_host_first => current
      xcfun_host_last => current
    ELSE
      xcfun_host_last%next => current
      xcfun_host_last => current
    ENDIF

    XCFUNfunctional = xc_new_functional()

    IF (XCFUNfunctional.EQ.-1) &
     &CALL LSQUIT('Error in xcfun_host_add_functional, xcfun internal maximum number of functionals reached',lupri)

    admm_ggacorr = .FALSE.
    IPOS = INDEX(DFTfuncString,'GGAKEY')
    if (IPOS.NE.0) then
       !GGAkeyString = .TRUE.
       call determine_nStrings(DFTfuncString(IPOS+7:),nStrings)
       
       allocate(DFTfuncStringSingle(nStrings))
       allocate(WeightSingle(nStrings))
       call trim_strings(DFTfuncString(IPOS+7:),nStrings,&
            & DFTfuncStringSingle,WeightSingle)
       
    else if ((INDEX(DFTfuncString,'LDAERF')).NE.0) then
       nStrings = 3
       !ErfString = .true.
       allocate(DFTfuncStringSingle(nStrings))
       allocate(WeightSingle(nStrings))
       call trim_strings(DFTfuncString,nStrings,&
            & DFTfuncStringSingle,WeightSingle)
       WeightSingle(1) = 1.0E0_realk
!!$    else if ((INDEX(DFTfuncString,'CAM')).NE.0) then
!!$       !FIX ME: ALPHA; BETA; MU
!!$       call determine_nStrings(DFTfuncString,nStrings)
!!$       allocate(DFTfuncStringSingle(nStrings))
!!$       allocate(WeightSingle(nStrings))
!!$       if (nStrings==1) then
!!$          call trim_strings(DFTfuncString,nStrings,&
!!$               & DFTfuncStringSingle,WeightSingle)
!!$          WeightSingle(1) = 1.0E0_realk
!!$       else
!!$          call lsquit('It is not possible to tune CAMB3LYP with XCFUN yet',lupri)
!!$       endif
    else
       call determine_nStrings(DFTfuncString,nStrings)
       
       allocate(DFTfuncStringSingle(nStrings))
       allocate(WeightSingle(nStrings))
       call trim_strings(DFTfuncString,nStrings,DFTfuncStringSingle,WeightSingle)
       do I=1,nStrings
          admm_ggacorr_single = ( (INDEX(DFTfuncStringSingle(I),'B88X')).NE.0 ) .OR.    &
         &                      ( (INDEX(DFTfuncStringSingle(I),'LDAX')).NE.0 ) .OR.    &
         &                      ( (INDEX(DFTfuncStringSingle(I),'PBEX')).NE.0 ) .OR.    &
         &                      ( (INDEX(DFTfuncStringSingle(I),'KT3X')).NE.0 ) .OR.    &
         &                      ( (INDEX(DFTfuncStringSingle(I),'OPTX')).NE.0 ) .OR.    &
         &                      ( (INDEX(DFTfuncStringSingle(I),'CAMCOMPX')).NE.0) .OR. &
         &                      ( (INDEX(DFTfuncStringSingle(I),'tfk')).NE.0)
          IF (admm_ggacorr_single) THEN
            WeightSingle(I) = hfweight
            admm_ggacorr = .TRUE.
          ELSE
            WeightSingle(I) = 1.0E0_realk
          ENDIF
       enddo
    ENDIF

    

    !XCFUNfunctional = xc_new_functional()
    do I=1,nStrings
       ierr = xc_set(XCFUNfunctional,DFTfuncStringSingle(I),WeightSingle(I))
       IF(ierr.NE.0)THEN
          print*,'The functional or parameter name:',DFTfuncStringSingle(I),'was not recognized'
          print*,' by the XCFUN program '
          write(lupri,*)'The functional or parameter name:',DFTfuncStringSingle(I),'was not recognized'
          write(lupri,*)' by the XCFUN program '
          call lsquit('xcfun_host_add_functional: error',lupri)
       ENDIF
    enddo   
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
    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
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
             XCFUN_DOGGA = .FALSE.
             XCFUN_DOMETA = .TRUE.
             call lsquit('Meta functionals current do not work.',lupri)
          ENDIF
       ELSE
          WRITE(lupri,*)'The Functional chosen is a GGA type functional'
          XCFUN_DOGGA = .TRUE.
          XCFUN_DOMETA = .FALSE.
       ENDIF
    ELSE
       WRITE(lupri,*)'The Functional chosen is a LDA type functional'
       XCFUN_DOGGA = .FALSE.
       XCFUN_DOMETA = .FALSE.
    ENDIF
    IF (admm_ggacorr) THEN
      !hfweight is an input, not an output
    ELSE
      ierrHF = xc_get(XCFUNfunctional,'EXX',hfweight)
      IF(ierrHF.NE.0)THEN
         call lsquit('Error determining exact (HF) exchange weight.',lupri)
      ENDIF
      IF(ABS(hfweight).GT.1.0E-16_realk)THEN
         WRITE(lupri,*)'The Functional chosen contains an exact exchange contribution'
         WRITE(lupri,*)'with the weight:', hfweight
      ELSE
         WRITE(lupri,*)'The Functional chosen contains no exact exchange contribution'
      ENDIF
    ENDIF
    deallocate(DFTfuncStringSingle)
    deallocate(WeightSingle)

    num_xc_functionals      = num_xc_functionals + 1
    current%DFTfuncString   = DFTfuncString
    current%XCFUNfunctional = XCFUNfunctional
    current%hfweight        = hfweight
    current%gga             = XCFUN_DOGGA
    current%meta            = XCFUN_DOMETA
    current%len             = len_trim(DFTfuncString)
    xcfun_host_current     => current
    nullify(current%next)

  end subroutine xcfun_host_add_functional
#endif

  subroutine xcfun_host_augment_func(DFTfuncString,hfweight,lupri)
    implicit none
    integer                     :: lupri
    character(*),intent(in)     :: DFTfuncString
    real(realk),intent(in)      :: hfweight
#ifdef VAR_XCFUN
    !
    type(xc_functional),pointer :: current
    character(len=1024)         :: augmentedFunc
    integer                     :: l1,l2
    real(realk)                 :: weight

    IF (.NOT.associated(xcfun_host_current)) CALL LSQUIT('Error in xcfun_host_augment_func. Current functional not set',-1)
    current => xcfun_host_current
    l1=current%len
    l2=LEN_TRIM(DFTfuncString)
    write(augmentedFunc,'(1024X)')
    augmentedFunc(1:l1) = current%DFTfuncString(1:l1)
    augmentedFunc(l1+1:l1+1) = ' '
    augmentedFunc(l1+2:l1+1+l2) = DFTfuncString(1:l2)

    weight = hfweight

    call xcfun_host_init(augmentedFunc,weight,lupri)

#endif
  end subroutine xcfun_host_augment_func


  subroutine xcfundftreport(lupri)
    implicit none
    integer :: lupri
#ifdef VAR_XCFUN
    integer     :: Ipos,nStrings,ierr,I,ierrLDA,ierrGGA,ierrMETA,ierrHF
    logical     :: GGAkeyString
    real(realk) :: hfweight
    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
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
             XCFUN_DOGGA = .FALSE.
             XCFUN_DOMETA = .TRUE.
             call lsquit('Meta functionals current do not work.',lupri)
          ENDIF
       ELSE
          WRITE(lupri,*)'The Functional chosen is a GGA type functional'
          XCFUN_DOGGA = .TRUE.
          XCFUN_DOMETA = .FALSE.
       ENDIF
    ELSE
       WRITE(lupri,*)'The Functional chosen is a LDA type functional'
       XCFUN_DOGGA = .FALSE.
       XCFUN_DOMETA = .FALSE.
    ENDIF
    ierrHF = xc_get(XCFUNfunctional,'EXX',hfweight)
    IF(ierrHF.NE.0)THEN
       call lsquit('Error determining exact (HF) exchange weight.',lupri)
    ENDIF
    IF(ABS(hfweight).GT.1.0E-16_realk)THEN
       WRITE(lupri,*)'The Functional chosen contains a exact exchange contribution'
       WRITE(lupri,*)'with the weight:', hfweight
    ELSE
       WRITE(lupri,*)'The Functional chosen contains no exact exchange contribution'
    ENDIF
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfundftreport

  subroutine xcfun_host_type(DOGGA,DOMETA)
    implicit none
    LOGICAL,intent(inout) :: DOGGA,DOMETA
    !
    integer :: Ipos,nStrings,ierr,I,ierrLDA,ierrGGA,ierrMETA,ierrHF
    logical :: GGAkeyString
    real(realk),pointer :: WeightSingle(:)
#ifdef VAR_XCFUN
    DOGGA = XCFUN_DOGGA
    DOMETA = XCFUN_DOMETA
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_host_type

  subroutine xcfun_host_set_order(order,unres,type)
    implicit none
    !> to which order the partial derivatives needs to be calculated
    integer,intent(in) :: order
    !> unrestricted or closed shell
    logical,intent(in) :: unres
    !> type : LDA,GGA,META see xcfun_type_lda,xcfun_type_gga,xcfun_type_metagga
    integer,intent(in) :: type
    integer :: ierrLDA,ierrGGA,ierrMETA
#ifdef VAR_XCFUN    
    IF(type.EQ.xcfun_type_lda)THEN
       if(unres)THEN
          ierrLDA = xc_eval_setup(XCFUNfunctional,XC_A_B,XC_PARTIAL_DERIVATIVES,order)
       else
          ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,order)
       endif
       IF(ierrLDA.NE.0) then
          print*,'ierrLDA',ierrLDA,'order',order,'unres',unres
          print*,'error in xc_eval_setup() for LDA'
          call lsquit('xcfun_set_order LDA error',-1)
       endif
    ELSEIF(type.EQ.xcfun_type_gga)THEN
       if(unres)THEN
          ierrGGA = xc_eval_setup(XCFUNfunctional,XC_A_B_GAA_GAB_GBB,XC_PARTIAL_DERIVATIVES,order)
       else
          ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ,XC_PARTIAL_DERIVATIVES,order)
       endif
       IF(ierrGGA.NE.0) then
          print*,'ierrGGA',ierrGGA,'order',order,'unres',unres
          print*,'error in xc_eval_setup() for GGA'
          call lsquit('xcfun_set_order GGA error',-1)
       endif
    ELSEIF(type.EQ.xcfun_type_metagga)THEN
       call lsquit('METAGGA not yet working for xcfun/lsdalton combination',-1)
       if(unres)THEN
!          ierrMETA = xc_eval_setup(XCFUNfunctional,XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB,XC_PARTIAL_DERIVATIVES,order)
       else
!          ierrMETA = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ_TAUN,XC_PARTIAL_DERIVATIVES,order)
       endif
       IF(ierrMETA.NE.0) then
          print*,'ierrMETA',ierrMETA,'order',order,'unres',unres
          print*,'error in xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,order)'
          call lsquit('xcfun_set_order METAGGA error',-1)
       endif
    ELSE
       call lsquit('unknown xcfun type in xcfun_set_order',-1)
    ENDIF
#endif
  end subroutine xcfun_host_set_order
  
!LDA PART

  subroutine xcfun_lda_xc_eval(XCFUNINPUT,XCFUNOUTPUT,NPNT)
    implicit none
    INTEGER,intent(IN)        :: NPNT
    REAL(REALK),intent(in)    :: XCFUNINPUT(1,NPNT)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(2,NPNT)
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    call xc_eval(XCFUNfunctional,NPNT,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_lda_xc_eval

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

  subroutine xcfun2_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(1,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(3,1)
    integer :: ierrLDA
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
!    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,2)
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
!    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d^2 Exc/d rho^2
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun2_lda_xc_single_eval

  subroutine xcfun3_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(1,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(4,1)
    integer :: ierrLDA
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
!    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,2)
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
!    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d^2 Exc/d rho^2
    ! XCFUNOUTPUT(4,1) - d^3 Exc/d rho^3
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun3_lda_xc_single_eval

!GGA PART

  subroutine xcfun_gga_unres_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT, ORDER)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(8,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(165,1)
    integer, intent(in) :: order
! Input:
    !rho_a   = XCFUNINPUT(1,1) 
    !rho_b   = XCFUNINPUT(2,1) 
    !agrad_x = XCFUNINPUT(3,1) 
    !agrad_y = XCFUNINPUT(4,1) 
    !agrad_z = XCFUNINPUT(5,1) 
    !bgrad_x = XCFUNINPUT(6,1) 
    !bgrad_y = XCFUNINPUT(7,1) 
    !bgrad_z = XCFUNINPUT(8,1)
! Output
    ! Order 0
    ! out(1,1) Exc
    ! Order 1
    ! out(2,1) d^1 Exc / d na
    ! out(3,1) d^1 Exc / d nb
    ! out(4,1) d^1 Exc / d ax
    ! out(5,1) d^1 Exc / d ay
    ! out(6,1) d^1 Exc / d az
    ! out(7,1) d^1 Exc / d bx
    ! out(8,1) d^1 Exc / d by
    ! out(9,1) d^1 Exc / d bz
    ! Order 2
    ! out(10,1) d^2 Exc / d na na
    ! out(11,1) d^2 Exc / d na nb
    ! out(12,1) d^2 Exc / d na ax
    ! out(13,1) d^2 Exc / d na ay
    ! out(14,1) d^2 Exc / d na az
    ! out(15,1) d^2 Exc / d na bx
    ! out(16,1) d^2 Exc / d na by
    ! out(17,1) d^2 Exc / d na bz
    ! out(18,1) d^2 Exc / d nb nb
    ! out(19,1) d^2 Exc / d nb ax
    ! out(20,1) d^2 Exc / d nb ay
    ! out(21,1) d^2 Exc / d nb az
    ! out(22,1) d^2 Exc / d nb bx
    ! out(23,1) d^2 Exc / d nb by
    ! out(24,1) d^2 Exc / d nb bz
    ! out(25,1) d^2 Exc / d ax ax
    ! out(26,1) d^2 Exc / d ax ay
    ! out(27,1) d^2 Exc / d ax az
    ! out(28,1) d^2 Exc / d ax bx
    ! out(29,1) d^2 Exc / d ax by
    ! out(30,1) d^2 Exc / d ax bz
    ! out(31,1) d^2 Exc / d ay ay
    ! out(32,1) d^2 Exc / d ay az
    ! out(33,1) d^2 Exc / d ay bx
    ! out(34,1) d^2 Exc / d ay by
    ! out(35,1) d^2 Exc / d ay bz
    ! out(36,1) d^2 Exc / d az az
    ! out(37,1) d^2 Exc / d az bx
    ! out(38,1) d^2 Exc / d az by
    ! out(39,1) d^2 Exc / d az bz
    ! out(40,1) d^2 Exc / d bx bx
    ! out(41,1) d^2 Exc / d bx by
    ! out(42,1) d^2 Exc / d bx bz
    ! out(43,1) d^2 Exc / d by by
    ! out(44,1) d^2 Exc / d by bz
    ! out(45,1) d^2 Exc / d bz bz
    ! Order 3
    ! out(46,1) d^3 Exc / d na na na
    ! out(47,1) d^3 Exc / d na na nb
    ! out(48,1) d^3 Exc / d na na ax
    ! out(49,1) d^3 Exc / d na na ay
    ! out(50,1) d^3 Exc / d na na az
    ! out(51,1) d^3 Exc / d na na bx
    ! out(52,1) d^3 Exc / d na na by
    ! out(53,1) d^3 Exc / d na na bz
    ! out(54,1) d^3 Exc / d na nb nb
    ! out(55,1) d^3 Exc / d na nb ax
    ! out(56,1) d^3 Exc / d na nb ay
    ! out(57,1) d^3 Exc / d na nb az
    ! out(58,1) d^3 Exc / d na nb bx
    ! out(59,1) d^3 Exc / d na nb by
    ! out(60,1) d^3 Exc / d na nb bz
    ! out(61,1) d^3 Exc / d na ax ax
    ! out(62,1) d^3 Exc / d na ax ay
    ! out(63,1) d^3 Exc / d na ax az
    ! out(64,1) d^3 Exc / d na ax bx
    ! out(65,1) d^3 Exc / d na ax by
    ! out(66,1) d^3 Exc / d na ax bz
    ! out(67,1) d^3 Exc / d na ay ay
    ! out(68,1) d^3 Exc / d na ay az
    ! out(69,1) d^3 Exc / d na ay bx
    ! out(70,1) d^3 Exc / d na ay by
    ! out(71,1) d^3 Exc / d na ay bz
    ! out(72,1) d^3 Exc / d na az az
    ! out(73,1) d^3 Exc / d na az bx
    ! out(74,1) d^3 Exc / d na az by
    ! out(75,1) d^3 Exc / d na az bz
    ! out(76,1) d^3 Exc / d na bx bx
    ! out(77,1) d^3 Exc / d na bx by
    ! out(78,1) d^3 Exc / d na bx bz
    ! out(79,1) d^3 Exc / d na by by
    ! out(80,1) d^3 Exc / d na by bz
    ! out(81,1) d^3 Exc / d na bz bz
    ! out(82,1) d^3 Exc / d nb nb nb
    ! out(83,1) d^3 Exc / d nb nb ax
    ! out(84,1) d^3 Exc / d nb nb ay
    ! out(85,1) d^3 Exc / d nb nb az
    ! out(86,1) d^3 Exc / d nb nb bx
    ! out(87,1) d^3 Exc / d nb nb by
    ! out(88,1) d^3 Exc / d nb nb bz
    ! out(89,1) d^3 Exc / d nb ax ax
    ! out(90,1) d^3 Exc / d nb ax ay
    ! out(91,1) d^3 Exc / d nb ax az
    ! out(92,1) d^3 Exc / d nb ax bx
    ! out(93,1) d^3 Exc / d nb ax by
    ! out(94,1) d^3 Exc / d nb ax bz
    ! out(95,1) d^3 Exc / d nb ay ay
    ! out(96,1) d^3 Exc / d nb ay az
    ! out(97,1) d^3 Exc / d nb ay bx
    ! out(98,1) d^3 Exc / d nb ay by
    ! out(99,1) d^3 Exc / d nb ay bz
    ! out(100,1) d^3 Exc / d nb az az
    ! out(101,1) d^3 Exc / d nb az bx
    ! out(102,1) d^3 Exc / d nb az by
    ! out(103,1) d^3 Exc / d nb az bz
    ! out(104,1) d^3 Exc / d nb bx bx
    ! out(105,1) d^3 Exc / d nb bx by
    ! out(106,1) d^3 Exc / d nb bx bz
    ! out(107,1) d^3 Exc / d nb by by
    ! out(108,1) d^3 Exc / d nb by bz
    ! out(109,1) d^3 Exc / d nb bz bz
    ! out(110,1) d^3 Exc / d ax ax ax
    ! out(111,1) d^3 Exc / d ax ax ay
    ! out(112,1) d^3 Exc / d ax ax az
    ! out(113,1) d^3 Exc / d ax ax bx
    ! out(114,1) d^3 Exc / d ax ax by
    ! out(115,1) d^3 Exc / d ax ax bz
    ! out(116,1) d^3 Exc / d ax ay ay
    ! out(117,1) d^3 Exc / d ax ay az
    ! out(118,1) d^3 Exc / d ax ay bx
    ! out(119,1) d^3 Exc / d ax ay by
    ! out(120,1) d^3 Exc / d ax ay bz
    ! out(121,1) d^3 Exc / d ax az az
    ! out(122,1) d^3 Exc / d ax az bx
    ! out(123,1) d^3 Exc / d ax az by
    ! out(124,1) d^3 Exc / d ax az bz
    ! out(125,1) d^3 Exc / d ax bx bx
    ! out(126,1) d^3 Exc / d ax bx by
    ! out(127,1) d^3 Exc / d ax bx bz
    ! out(128,1) d^3 Exc / d ax by by
    ! out(129,1) d^3 Exc / d ax by bz
    ! out(130,1) d^3 Exc / d ax bz bz
    ! out(131,1) d^3 Exc / d ay ay ay
    ! out(132,1) d^3 Exc / d ay ay az
    ! out(133,1) d^3 Exc / d ay ay bx
    ! out(134,1) d^3 Exc / d ay ay by
    ! out(135,1) d^3 Exc / d ay ay bz
    ! out(136,1) d^3 Exc / d ay az az
    ! out(137,1) d^3 Exc / d ay az bx
    ! out(138,1) d^3 Exc / d ay az by
    ! out(139,1) d^3 Exc / d ay az bz
    ! out(140,1) d^3 Exc / d ay bx bx
    ! out(141,1) d^3 Exc / d ay bx by
    ! out(142,1) d^3 Exc / d ay bx bz
    ! out(143,1) d^3 Exc / d ay by by
    ! out(144,1) d^3 Exc / d ay by bz
    ! out(145,1) d^3 Exc / d ay bz bz
    ! out(146,1) d^3 Exc / d az az az
    ! out(147,1) d^3 Exc / d az az bx
    ! out(148,1) d^3 Exc / d az az by
    ! out(149,1) d^3 Exc / d az az bz
    ! out(150,1) d^3 Exc / d az bx bx
    ! out(151,1) d^3 Exc / d az bx by
    ! out(152,1) d^3 Exc / d az bx bz
    ! out(153,1) d^3 Exc / d az by by
    ! out(154,1) d^3 Exc / d az by bz
    ! out(155,1) d^3 Exc / d az bz bz
    ! out(156,1) d^3 Exc / d bx bx bx
    ! out(157,1) d^3 Exc / d bx bx by
    ! out(158,1) d^3 Exc / d bx bx bz
    ! out(159,1) d^3 Exc / d bx by by
    ! out(160,1) d^3 Exc / d bx by bz
    ! out(161,1) d^3 Exc / d bx bz bz
    ! out(162,1) d^3 Exc / d by by by
    ! out(163,1) d^3 Exc / d by by bz
    ! out(164,1) d^3 Exc / d by bz bz
    ! out(165,1) d^3 Exc / d bz bz bz
    integer :: ierr
#ifdef VAR_XCFUN
    ierr = xc_eval_setup(XCFUNfunctional,XC_A_B_AX_AY_AZ_BX_BY_BZ,XC_PARTIAL_DERIVATIVES,order)
    IF(ierr.NE.0) then
       print*,'ierr from xcfun',ierr
       call lsquit('Unexpected error from xcfun',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_gga_unres_components_xc_single_eval

  subroutine xcfun_gga_components_xc_single_eval(XCFUNINPUT,ndim,XCFUNOUTPUT,ORDER)
    implicit none
    integer, intent(in) :: order,ndim
!ndim depend on order 
!order 0 => dim = 1
!order 1 => dim = 5
!order 2 => dim = 35
    REAL(REALK),intent(in) :: XCFUNINPUT(4,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(ndim,1)
! Input:
    !rho   = XCFUNINPUT(1,1)
    !grad_x = XCFUNINPUT(2,1)
    !grad_y = XCFUNINPUT(3,1)
    !grad_z = XCFUNINPUT(4,1)
! Output
    ! Order 0
    ! out(1,1) Exc
    ! Order 1
    ! out(2,1) d^1 Exc / d n
    ! out(3,1) d^1 Exc / d nx
    ! out(4,1) d^1 Exc / d ny
    ! out(5,1) d^1 Exc / d nz
    ! Order 2
    ! out(6,1) d^2 Exc / d n n
    ! out(7,1) d^2 Exc / d n nx
    ! out(8,1) d^2 Exc / d n ny
    ! out(9,1) d^2 Exc / d n nz
    ! out(10,1) d^2 Exc / d nx nx
    ! out(11,1) d^2 Exc / d nx ny
    ! out(12,1) d^2 Exc / d nx nz
    ! out(13,1) d^2 Exc / d ny ny
    ! out(14,1) d^2 Exc / d ny nz
    ! out(15,1) d^2 Exc / d nz nz
    ! Order 3
    ! out(16,1) d^3 Exc / d n n n
    ! out(17,1) d^3 Exc / d n n nx
    ! out(18,1) d^3 Exc / d n n ny
    ! out(19,1) d^3 Exc / d n n nz
    ! out(20,1) d^3 Exc / d n nx nx
    ! out(21,1) d^3 Exc / d n nx ny
    ! out(22,1) d^3 Exc / d n nx nz
    ! out(23,1) d^3 Exc / d n ny ny
    ! out(24,1) d^3 Exc / d n ny nz
    ! out(25,1) d^3 Exc / d n nz nz
    ! out(26,1) d^3 Exc / d nx nx nx
    ! out(27,1) d^3 Exc / d nx nx ny
    ! out(28,1) d^3 Exc / d nx nx nz
    ! out(29,1) d^3 Exc / d nx ny ny
    ! out(30,1) d^3 Exc / d nx ny nz
    ! out(31,1) d^3 Exc / d nx nz nz
    ! out(32,1) d^3 Exc / d ny ny ny
    ! out(33,1) d^3 Exc / d ny ny nz
    ! out(34,1) d^3 Exc / d ny nz nz
    ! out(35,1) d^3 Exc / d nz nz nz
    integer :: ierr
#ifdef VAR_XCFUN
    ierr = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ,XC_PARTIAL_DERIVATIVES,order)
    IF(ierr.NE.0) then
       print*,'ierr from xcfun',ierr
       call lsquit('Unexpected error from xcfun',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_gga_components_xc_single_eval


  subroutine xcfun_gga_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(3,1)
#ifdef VAR_XCFUN
    integer :: ierrGGA
    !rho = XCFUNINPUT(1,1) 
    !gg = |grad|^2 = XCFUNINPUT(2,1)
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,1)
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ! XCFUNOUTPUT(1,1) - Exc
    ! XCFUNOUTPUT(2,1) - d Exc/d rho
    ! XCFUNOUTPUT(3,1) - d Exc/d gg
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_gga_xc_single_eval

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

!METAGGA






  subroutine xcfun_metagga_unres_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT, ORDER)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(10,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(286,1)
    integer, intent(in) :: order
! Input:
    !rho_a   = XCFUNINPUT(1,1) 
    !rho_b   = XCFUNINPUT(2,1) 
    !agrad_x = XCFUNINPUT(3,1) 
    !agrad_y = XCFUNINPUT(4,1) 
    !agrad_z = XCFUNINPUT(5,1) 
    !bgrad_x = XCFUNINPUT(6,1) 
    !bgrad_y = XCFUNINPUT(7,1) 
    !bgrad_z = XCFUNINPUT(8,1) 
    !taua    = XCFUNINPUT(9,1) 
    !taub    = XCFUNINPUT(10,1) 
! Output
    ! Order 0
    ! out(1,1) Exc
    ! Order 1
    ! out(2,1) d^1 Exc / d na
    ! out(3,1) d^1 Exc / d nb
    ! out(4,1) d^1 Exc / d ax
    ! out(5,1) d^1 Exc / d ay
    ! out(6,1) d^1 Exc / d az
    ! out(7,1) d^1 Exc / d bx
    ! out(8,1) d^1 Exc / d by
    ! out(9,1) d^1 Exc / d bz
    ! out(10,1) d^1 Exc / d taua
    ! out(11,1) d^1 Exc / d taub
    ! Order 2
    ! out(12,1) d^2 Exc / d na na
    ! out(13,1) d^2 Exc / d na nb
    ! out(14,1) d^2 Exc / d na ax
    ! out(15,1) d^2 Exc / d na ay
    ! out(16,1) d^2 Exc / d na az
    ! out(17,1) d^2 Exc / d na bx
    ! out(18,1) d^2 Exc / d na by
    ! out(19,1) d^2 Exc / d na bz
    ! out(20,1) d^2 Exc / d na taua
    ! out(21,1) d^2 Exc / d na taub
    ! out(22,1) d^2 Exc / d nb nb
    ! out(23,1) d^2 Exc / d nb ax
    ! out(24,1) d^2 Exc / d nb ay
    ! out(25,1) d^2 Exc / d nb az
    ! out(26,1) d^2 Exc / d nb bx
    ! out(27,1) d^2 Exc / d nb by
    ! out(28,1) d^2 Exc / d nb bz
    ! out(29,1) d^2 Exc / d nb taua
    ! out(30,1) d^2 Exc / d nb taub
    ! out(31,1) d^2 Exc / d ax ax
    ! out(32,1) d^2 Exc / d ax ay
    ! out(33,1) d^2 Exc / d ax az
    ! out(34,1) d^2 Exc / d ax bx
    ! out(35,1) d^2 Exc / d ax by
    ! out(36,1) d^2 Exc / d ax bz
    ! out(37,1) d^2 Exc / d ax taua
    ! out(38,1) d^2 Exc / d ax taub
    ! out(39,1) d^2 Exc / d ay ay
    ! out(40,1) d^2 Exc / d ay az
    ! out(41,1) d^2 Exc / d ay bx
    ! out(42,1) d^2 Exc / d ay by
    ! out(43,1) d^2 Exc / d ay bz
    ! out(44,1) d^2 Exc / d ay taua
    ! out(45,1) d^2 Exc / d ay taub
    ! out(46,1) d^2 Exc / d az az
    ! out(47,1) d^2 Exc / d az bx
    ! out(48,1) d^2 Exc / d az by
    ! out(49,1) d^2 Exc / d az bz
    ! out(50,1) d^2 Exc / d az taua
    ! out(51,1) d^2 Exc / d az taub
    ! out(52,1) d^2 Exc / d bx bx
    ! out(53,1) d^2 Exc / d bx by
    ! out(54,1) d^2 Exc / d bx bz
    ! out(55,1) d^2 Exc / d bx taua
    ! out(56,1) d^2 Exc / d bx taub
    ! out(57,1) d^2 Exc / d by by
    ! out(58,1) d^2 Exc / d by bz
    ! out(59,1) d^2 Exc / d by taua
    ! out(60,1) d^2 Exc / d by taub
    ! out(61,1) d^2 Exc / d bz bz
    ! out(62,1) d^2 Exc / d bz taua
    ! out(63,1) d^2 Exc / d bz taub
    ! out(64,1) d^2 Exc / d taua taua
    ! out(65,1) d^2 Exc / d taua taub
    ! out(66,1) d^2 Exc / d taub taub
    ! Order 3
    ! out(67,1) d^3 Exc / d na na na
    ! out(68,1) d^3 Exc / d na na nb
    ! out(69,1) d^3 Exc / d na na ax
    ! out(70,1) d^3 Exc / d na na ay
    ! out(71,1) d^3 Exc / d na na az
    ! out(72,1) d^3 Exc / d na na bx
    ! out(73,1) d^3 Exc / d na na by
    ! out(74,1) d^3 Exc / d na na bz
    ! out(75,1) d^3 Exc / d na na taua
    ! out(76,1) d^3 Exc / d na na taub
    ! out(77,1) d^3 Exc / d na nb nb
    ! out(78,1) d^3 Exc / d na nb ax
    ! out(79,1) d^3 Exc / d na nb ay
    ! out(80,1) d^3 Exc / d na nb az
    ! out(81,1) d^3 Exc / d na nb bx
    ! out(82,1) d^3 Exc / d na nb by
    ! out(83,1) d^3 Exc / d na nb bz
    ! out(84,1) d^3 Exc / d na nb taua
    ! out(85,1) d^3 Exc / d na nb taub
    ! out(86,1) d^3 Exc / d na ax ax
    ! out(87,1) d^3 Exc / d na ax ay
    ! out(88,1) d^3 Exc / d na ax az
    ! out(89,1) d^3 Exc / d na ax bx
    ! out(90,1) d^3 Exc / d na ax by
    ! out(91,1) d^3 Exc / d na ax bz
    ! out(92,1) d^3 Exc / d na ax taua
    ! out(93,1) d^3 Exc / d na ax taub
    ! out(94,1) d^3 Exc / d na ay ay
    ! out(95,1) d^3 Exc / d na ay az
    ! out(96,1) d^3 Exc / d na ay bx
    ! out(97,1) d^3 Exc / d na ay by
    ! out(98,1) d^3 Exc / d na ay bz
    ! out(99,1) d^3 Exc / d na ay taua
    ! out(100,1) d^3 Exc / d na ay taub
    ! out(101,1) d^3 Exc / d na az az
    ! out(102,1) d^3 Exc / d na az bx
    ! out(103,1) d^3 Exc / d na az by
    ! out(104,1) d^3 Exc / d na az bz
    ! out(105,1) d^3 Exc / d na az taua
    ! out(106,1) d^3 Exc / d na az taub
    ! out(107,1) d^3 Exc / d na bx bx
    ! out(108,1) d^3 Exc / d na bx by
    ! out(109,1) d^3 Exc / d na bx bz
    ! out(110,1) d^3 Exc / d na bx taua
    ! out(111,1) d^3 Exc / d na bx taub
    ! out(112,1) d^3 Exc / d na by by
    ! out(113,1) d^3 Exc / d na by bz
    ! out(114,1) d^3 Exc / d na by taua
    ! out(115,1) d^3 Exc / d na by taub
    ! out(116,1) d^3 Exc / d na bz bz
    ! out(117,1) d^3 Exc / d na bz taua
    ! out(118,1) d^3 Exc / d na bz taub
    ! out(119,1) d^3 Exc / d na taua taua
    ! out(120,1) d^3 Exc / d na taua taub
    ! out(121,1) d^3 Exc / d na taub taub
    ! out(122,1) d^3 Exc / d nb nb nb
    ! out(123,1) d^3 Exc / d nb nb ax
    ! out(124,1) d^3 Exc / d nb nb ay
    ! out(125,1) d^3 Exc / d nb nb az
    ! out(126,1) d^3 Exc / d nb nb bx
    ! out(127,1) d^3 Exc / d nb nb by
    ! out(128,1) d^3 Exc / d nb nb bz
    ! out(129,1) d^3 Exc / d nb nb taua
    ! out(130,1) d^3 Exc / d nb nb taub
    ! out(131,1) d^3 Exc / d nb ax ax
    ! out(132,1) d^3 Exc / d nb ax ay
    ! out(133,1) d^3 Exc / d nb ax az
    ! out(134,1) d^3 Exc / d nb ax bx
    ! out(135,1) d^3 Exc / d nb ax by
    ! out(136,1) d^3 Exc / d nb ax bz
    ! out(137,1) d^3 Exc / d nb ax taua
    ! out(138,1) d^3 Exc / d nb ax taub
    ! out(139,1) d^3 Exc / d nb ay ay
    ! out(140,1) d^3 Exc / d nb ay az
    ! out(141,1) d^3 Exc / d nb ay bx
    ! out(142,1) d^3 Exc / d nb ay by
    ! out(143,1) d^3 Exc / d nb ay bz
    ! out(144,1) d^3 Exc / d nb ay taua
    ! out(145,1) d^3 Exc / d nb ay taub
    ! out(146,1) d^3 Exc / d nb az az
    ! out(147,1) d^3 Exc / d nb az bx
    ! out(148,1) d^3 Exc / d nb az by
    ! out(149,1) d^3 Exc / d nb az bz
    ! out(150,1) d^3 Exc / d nb az taua
    ! out(151,1) d^3 Exc / d nb az taub
    ! out(152,1) d^3 Exc / d nb bx bx
    ! out(153,1) d^3 Exc / d nb bx by
    ! out(154,1) d^3 Exc / d nb bx bz
    ! out(155,1) d^3 Exc / d nb bx taua
    ! out(156,1) d^3 Exc / d nb bx taub
    ! out(157,1) d^3 Exc / d nb by by
    ! out(158,1) d^3 Exc / d nb by bz
    ! out(159,1) d^3 Exc / d nb by taua
    ! out(160,1) d^3 Exc / d nb by taub
    ! out(161,1) d^3 Exc / d nb bz bz
    ! out(162,1) d^3 Exc / d nb bz taua
    ! out(163,1) d^3 Exc / d nb bz taub
    ! out(164,1) d^3 Exc / d nb taua taua
    ! out(165,1) d^3 Exc / d nb taua taub
    ! out(166,1) d^3 Exc / d nb taub taub
    ! out(167,1) d^3 Exc / d ax ax ax
    ! out(168,1) d^3 Exc / d ax ax ay
    ! out(169,1) d^3 Exc / d ax ax az
    ! out(170,1) d^3 Exc / d ax ax bx
    ! out(171,1) d^3 Exc / d ax ax by
    ! out(172,1) d^3 Exc / d ax ax bz
    ! out(173,1) d^3 Exc / d ax ax taua
    ! out(174,1) d^3 Exc / d ax ax taub
    ! out(175,1) d^3 Exc / d ax ay ay
    ! out(176,1) d^3 Exc / d ax ay az
    ! out(177,1) d^3 Exc / d ax ay bx
    ! out(178,1) d^3 Exc / d ax ay by
    ! out(179,1) d^3 Exc / d ax ay bz
    ! out(180,1) d^3 Exc / d ax ay taua
    ! out(181,1) d^3 Exc / d ax ay taub
    ! out(182,1) d^3 Exc / d ax az az
    ! out(183,1) d^3 Exc / d ax az bx
    ! out(184,1) d^3 Exc / d ax az by
    ! out(185,1) d^3 Exc / d ax az bz
    ! out(186,1) d^3 Exc / d ax az taua
    ! out(187,1) d^3 Exc / d ax az taub
    ! out(188,1) d^3 Exc / d ax bx bx
    ! out(189,1) d^3 Exc / d ax bx by
    ! out(190,1) d^3 Exc / d ax bx bz
    ! out(191,1) d^3 Exc / d ax bx taua
    ! out(192,1) d^3 Exc / d ax bx taub
    ! out(193,1) d^3 Exc / d ax by by
    ! out(194,1) d^3 Exc / d ax by bz
    ! out(195,1) d^3 Exc / d ax by taua
    ! out(196,1) d^3 Exc / d ax by taub
    ! out(197,1) d^3 Exc / d ax bz bz
    ! out(198,1) d^3 Exc / d ax bz taua
    ! out(199,1) d^3 Exc / d ax bz taub
    ! out(200,1) d^3 Exc / d ax taua taua
    ! out(201,1) d^3 Exc / d ax taua taub
    ! out(202,1) d^3 Exc / d ax taub taub
    ! out(203,1) d^3 Exc / d ay ay ay
    ! out(204,1) d^3 Exc / d ay ay az
    ! out(205,1) d^3 Exc / d ay ay bx
    ! out(206,1) d^3 Exc / d ay ay by
    ! out(207,1) d^3 Exc / d ay ay bz
    ! out(208,1) d^3 Exc / d ay ay taua
    ! out(209,1) d^3 Exc / d ay ay taub
    ! out(210,1) d^3 Exc / d ay az az
    ! out(211,1) d^3 Exc / d ay az bx
    ! out(212,1) d^3 Exc / d ay az by
    ! out(213,1) d^3 Exc / d ay az bz
    ! out(214,1) d^3 Exc / d ay az taua
    ! out(215,1) d^3 Exc / d ay az taub
    ! out(216,1) d^3 Exc / d ay bx bx
    ! out(217,1) d^3 Exc / d ay bx by
    ! out(218,1) d^3 Exc / d ay bx bz
    ! out(219,1) d^3 Exc / d ay bx taua
    ! out(220,1) d^3 Exc / d ay bx taub
    ! out(221,1) d^3 Exc / d ay by by
    ! out(222,1) d^3 Exc / d ay by bz
    ! out(223,1) d^3 Exc / d ay by taua
    ! out(224,1) d^3 Exc / d ay by taub
    ! out(225,1) d^3 Exc / d ay bz bz
    ! out(226,1) d^3 Exc / d ay bz taua
    ! out(227,1) d^3 Exc / d ay bz taub
    ! out(228,1) d^3 Exc / d ay taua taua
    ! out(229,1) d^3 Exc / d ay taua taub
    ! out(230,1) d^3 Exc / d ay taub taub
    ! out(231,1) d^3 Exc / d az az az
    ! out(232,1) d^3 Exc / d az az bx
    ! out(233,1) d^3 Exc / d az az by
    ! out(234,1) d^3 Exc / d az az bz
    ! out(235,1) d^3 Exc / d az az taua
    ! out(236,1) d^3 Exc / d az az taub
    ! out(237,1) d^3 Exc / d az bx bx
    ! out(238,1) d^3 Exc / d az bx by
    ! out(239,1) d^3 Exc / d az bx bz
    ! out(240,1) d^3 Exc / d az bx taua
    ! out(241,1) d^3 Exc / d az bx taub
    ! out(242,1) d^3 Exc / d az by by
    ! out(243,1) d^3 Exc / d az by bz
    ! out(244,1) d^3 Exc / d az by taua
    ! out(245,1) d^3 Exc / d az by taub
    ! out(246,1) d^3 Exc / d az bz bz
    ! out(247,1) d^3 Exc / d az bz taua
    ! out(248,1) d^3 Exc / d az bz taub
    ! out(249,1) d^3 Exc / d az taua taua
    ! out(250,1) d^3 Exc / d az taua taub
    ! out(251,1) d^3 Exc / d az taub taub
    ! out(252,1) d^3 Exc / d bx bx bx
    ! out(253,1) d^3 Exc / d bx bx by
    ! out(254,1) d^3 Exc / d bx bx bz
    ! out(255,1) d^3 Exc / d bx bx taua
    ! out(256,1) d^3 Exc / d bx bx taub
    ! out(257,1) d^3 Exc / d bx by by
    ! out(258,1) d^3 Exc / d bx by bz
    ! out(259,1) d^3 Exc / d bx by taua
    ! out(260,1) d^3 Exc / d bx by taub
    ! out(261,1) d^3 Exc / d bx bz bz
    ! out(262,1) d^3 Exc / d bx bz taua
    ! out(263,1) d^3 Exc / d bx bz taub
    ! out(264,1) d^3 Exc / d bx taua taua
    ! out(265,1) d^3 Exc / d bx taua taub
    ! out(266,1) d^3 Exc / d bx taub taub
    ! out(267,1) d^3 Exc / d by by by
    ! out(268,1) d^3 Exc / d by by bz
    ! out(269,1) d^3 Exc / d by by taua
    ! out(270,1) d^3 Exc / d by by taub
    ! out(271,1) d^3 Exc / d by bz bz
    ! out(272,1) d^3 Exc / d by bz taua
    ! out(273,1) d^3 Exc / d by bz taub
    ! out(274,1) d^3 Exc / d by taua taua
    ! out(275,1) d^3 Exc / d by taua taub
    ! out(276,1) d^3 Exc / d by taub taub
    ! out(277,1) d^3 Exc / d bz bz bz
    ! out(278,1) d^3 Exc / d bz bz taua
    ! out(279,1) d^3 Exc / d bz bz taub
    ! out(280,1) d^3 Exc / d bz taua taua
    ! out(281,1) d^3 Exc / d bz taua taub
    ! out(282,1) d^3 Exc / d bz taub taub
    ! out(283,1) d^3 Exc / d taua taua taua
    ! out(284,1) d^3 Exc / d taua taua taub
    ! out(285,1) d^3 Exc / d taua taub taub
    ! out(286,1) d^3 Exc / d taub taub taub
    integer :: ierr
#ifdef VAR_XCFUN
!    ierr = xc_eval_setup(XCFUNfunctional,XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB,XC_PARTIAL_DERIVATIVES,order)
!    IF(ierr.NE.0) then
!       print*,'ierr from xcfun',ierr
!       call lsquit('Unexpected error from xcfun',-1)
!    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_metagga_unres_components_xc_single_eval

  subroutine xcfun_metagga_components_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT,order)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(5,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(56,1)
    integer, intent(in) :: order
! Input
    !rho    = XCFUNINPUT(1,1) 
    !grad_x = XCFUNINPUT(2,1) 
    !grad_y = XCFUNINPUT(3,1) 
    !grad_z = XCFUNINPUT(4,1) 
    !tau    = XCFUNINPUT(5,1)
! Output
    ! Order 0
    ! out(1,1) Exc
    ! Order 1
    ! out(2,1) d^1 Exc / d n
    ! out(3,1) d^1 Exc / d nx
    ! out(4,1) d^1 Exc / d ny
    ! out(5,1) d^1 Exc / d nz
    ! out(6,1) d^1 Exc / d tau
    ! Order 2
    ! out(7,1) d^2 Exc / d n n
    ! out(8,1) d^2 Exc / d n nx
    ! out(9,1) d^2 Exc / d n ny
    ! out(10,1) d^2 Exc / d n nz
    ! out(11,1) d^2 Exc / d n tau
    ! out(12,1) d^2 Exc / d nx nx
    ! out(13,1) d^2 Exc / d nx ny
    ! out(14,1) d^2 Exc / d nx nz
    ! out(15,1) d^2 Exc / d nx tau
    ! out(16,1) d^2 Exc / d ny ny
    ! out(17,1) d^2 Exc / d ny nz
    ! out(18,1) d^2 Exc / d ny tau
    ! out(19,1) d^2 Exc / d nz nz
    ! out(20,1) d^2 Exc / d nz tau
    ! out(21,1) d^2 Exc / d tau tau
    ! Order 3
    ! out(22,1) d^3 Exc / d n n n
    ! out(23,1) d^3 Exc / d n n nx
    ! out(24,1) d^3 Exc / d n n ny
    ! out(25,1) d^3 Exc / d n n nz
    ! out(26,1) d^3 Exc / d n n tau
    ! out(27,1) d^3 Exc / d n nx nx
    ! out(28,1) d^3 Exc / d n nx ny
    ! out(29,1) d^3 Exc / d n nx nz
    ! out(30,1) d^3 Exc / d n nx tau
    ! out(31,1) d^3 Exc / d n ny ny
    ! out(32,1) d^3 Exc / d n ny nz
    ! out(33,1) d^3 Exc / d n ny tau
    ! out(34,1) d^3 Exc / d n nz nz
    ! out(35,1) d^3 Exc / d n nz tau
    ! out(36,1) d^3 Exc / d n tau tau
    ! out(37,1) d^3 Exc / d nx nx nx
    ! out(38,1) d^3 Exc / d nx nx ny
    ! out(39,1) d^3 Exc / d nx nx nz
    ! out(40,1) d^3 Exc / d nx nx tau
    ! out(41,1) d^3 Exc / d nx ny ny
    ! out(42,1) d^3 Exc / d nx ny nz
    ! out(43,1) d^3 Exc / d nx ny tau
    ! out(44,1) d^3 Exc / d nx nz nz
    ! out(45,1) d^3 Exc / d nx nz tau
    ! out(46,1) d^3 Exc / d nx tau tau
    ! out(47,1) d^3 Exc / d ny ny ny
    ! out(48,1) d^3 Exc / d ny ny nz
    ! out(49,1) d^3 Exc / d ny ny tau
    ! out(50,1) d^3 Exc / d ny nz nz
    ! out(51,1) d^3 Exc / d ny nz tau
    ! out(52,1) d^3 Exc / d ny tau tau
    ! out(53,1) d^3 Exc / d nz nz nz
    ! out(54,1) d^3 Exc / d nz nz tau
    ! out(55,1) d^3 Exc / d nz tau tau
    ! out(56,1) d^3 Exc / d tau tau tau
    integer :: ierr
#ifdef VAR_XCFUN
!    ierr = xc_eval_setup(XCFUNfunctional,XC_N_NX_NY_NZ_TAUN,XC_PARTIAL_DERIVATIVES,order)
!    IF(ierr.NE.0) then
!       print*,'ierr from xcfun',ierr
!       call lsquit('Unexpected error from xcfun',-1)
!    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_metagga_components_xc_single_eval

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

  subroutine xcfun_host_free()
    implicit none
#ifdef VAR_XCFUN
    call xc_free_functional(XCFUNfunctional)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_host_free

#ifdef VAR_XCFUN
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
                TMPweight(1:nweight) = string(i-nweight:i-1)
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

#endif

end module xcfun_host
