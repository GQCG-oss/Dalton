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

  subroutine xcfun_gga_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(3,1)
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    !|grad|^2 = XCFUNINPUT(2,1) 
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
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
    !|grad|^2 = XCFUNINPUT(2,1) = \sum_{\mu \nu} D_{\mu \nu} (\nabla chi_{mu} \times chi_{nu} + chi_{mu} \times \nabla chi_{nu})
    !tau = XCFUNINPUT(3,1)  = \sum_{\mu \nu} D_{\mu \nu} half nabla chi_{mu} \times \nabla chi_{nu}
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    !last is derivative wrt tau
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
    !|grad|^2 = XCFUNINPUT(2,1) 
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,2)
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,1)
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
    !|grad|^2 = XCFUNINPUT(2,1) 
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,3)
    IF(ierrGGA.NE.0) then
       print*,'ierrGGA',ierrGGA
       call lsquit('fun3 too large1',-1)
    endif
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,1)
    IF(ierrGGA.NE.0) call lsquit('fun3 too large2',-1)
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
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_lda_xc_single_eval

  subroutine xcfun2_lda_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(1,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(4,1)
    !
    integer :: ierrLDA
#ifdef VAR_XCFUN
    !rho = XCFUNINPUT(1,1) 
    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,2)
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun2_lda_xc_single_eval

  subroutine xcfun_gga_unres_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(5,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(6,1)
#ifdef VAR_XCFUN
    !rho_alpha = XCFUNINPUT(1,1) 
    !rho_beta = XCFUNINPUT(2,1) 
    !|grad_alpha|^2 = XCFUNINPUT(3,1) 
    !grad_alpha*grad_beta = XCFUNINPUT(4,1) 
    !|grad_beta|^2 = XCFUNINPUT(5,1) 
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)    
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
  end subroutine xcfun_gga_unres_xc_single_eval

  subroutine xcfun_lda_unres_xc_single_eval(XCFUNINPUT,XCFUNOUTPUT)
    implicit none
    REAL(REALK),intent(in) :: XCFUNINPUT(2,1)
    REAL(REALK),intent(inout) :: XCFUNOUTPUT(3,1)
#ifdef VAR_XCFUN
    !rho_alpha = XCFUNINPUT(1,1) 
    !rho_beta = XCFUNINPUT(2,1) 
    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
#else
    call lsquit('xcfun not activated -DVAR_XCFUN (can only be done using cmake)',-1)
#endif
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
