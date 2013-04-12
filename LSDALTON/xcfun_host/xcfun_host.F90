module xcfun_host
  use precision
!  use xcfun
  implicit none
  integer :: XCFunFunctional
  logical :: USEXCFUN
  public :: xcfun_host_init, xcfun_host_free, xcfun_gga_xc_single_eval,USEXCFUN
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
!    XCFUNfunctional = xc_new_functional()
    do I=1,nStrings
!       ierr = xc_set(XCFUNfunctional,DFTfuncStringSingle(I),WeightSingle(I))
       IF(ierr.NE.0)THEN
          print*,'The functional name:',DFTfuncStringSingle(I),'was not recognized'
          print*,' by the XCFUN program '
          write(lupri,*)'The functional name:',DFTfuncStringSingle(I),'was not recognized'
          write(lupri,*)' by the XCFUN program '
          call lsquit('xcfun_host_init: error',lupri)
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
!    ierrLDA = xc_eval_setup(XCFUNfunctional,XC_N,XC_PARTIAL_DERIVATIVES,1)
    ierrLDA=2
    IF(ierrLDA.NE.0)THEN
       !Test if this is a GGA type
!       ierrGGA = xc_eval_setup(XCFUNfunctional,XC_N_GNN,XC_PARTIAL_DERIVATIVES,1)
       IF(ierrGGA.NE.0)THEN
       !Test if this is a Meta type
!          ierrMETA = xc_eval_setup(XCFUNfunctional,XC_N_GNN_TAUN,XC_PARTIAL_DERIVATIVES,1)
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
!    ierrHF = xc_get(XCFUNfunctional,'EXX',hfweight)
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
!    call xc_eval(XCFUNfunctional,1,XCFUNINPUT,XCFUNOUTPUT)
  end subroutine xcfun_gga_xc_single_eval

  subroutine xcfun_host_free()
    implicit none
!    call xc_free_functional(XCFUNfunctional)
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
