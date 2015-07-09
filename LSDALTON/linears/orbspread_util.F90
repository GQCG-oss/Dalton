module orbspread_utilMod
!##########################################################
!#              ORBSPREAD MODULE                          #
!# Routines that are specific for orbspread localization  #
!# Routines called by solver (lin.trans. and precond.)    #
!# are not included in module.                            #
!#                                                        #
!##########################################################
use precision
use loc_utils
use typedef
use typedeftype
use loc_types
use matrix_module, only: matrix
use matrix_operations 
use matrix_operations_aux
use LSTIMING
use integralInterfaceMod
use decompMod !orbspread_data
use arhDensity

contains
  subroutine orbspread_free(orbspread_input,DoNotAllocateTmpM)
    implicit none
    type(orbspread_data) :: orbspread_input
    logical,optional :: DoNotAllocateTmpM
    integer :: i
    logical :: DeAllocTMP

    DeAllocTMP = .TRUE.
    IF(present(DoNotAllocateTmpM))THEN
       DeAllocTMP = .NOT.DoNotAllocateTmpM
    ENDIF


    call mat_free(orbspread_input%Q)

    do i=1,3
       call mat_free(orbspread_input%R(i))
    enddo
    IF(DeAllocTMP)THEN
       call mem_dealloc(orbspread_input%spread2) 
!       do i=1,4
!          call mat_free(orbspread_input%tmpM(i))
!       enddo
    ENDIF

  end subroutine orbspread_free

  subroutine orbspread_propint_free(orbspread_input)
    implicit none
    type(orbspread_data) :: orbspread_input
    integer :: i

    do i=2,5
       call mat_free(orbspread_input%propint(i))
    enddo

  end subroutine orbspread_propint_free


  subroutine orbspread_propint(orbspread_input,ls)
    implicit none
    type(orbspread_data), intent(inout) :: orbspread_input
    TYPE(lsitem), intent(inout) :: ls !ls%setting is changed in the II_something routines

    integer, parameter  :: nderiv=2, nMAT=10
    integer i,n

    ! init and compute propint property integrals
    ! propint(2:4) -> DIPX,DIPY,DIPZ
    ! propint(5)   -> SECX+SECY+SECZ 

    !     n   = ls%setting%BASIS(1)%p%REGULAR%nbast
    n = getNbasis(AORdefault,Contractedinttype,ls%SETTING%MOLECULE(1)%p,ls%lupri)

    IF(.FALSE.)THEN !speed optimized version
       DO i=1,nMAT
          call mat_init(orbspread_input%propint(I),n,n)
       ENDDO
       call II_get_carmom(6,6,ls%setting,orbspread_input%propint,nMAT,nderiv,0E0_realk,0E0_realk,0E0_realk)
       call mat_free(orbspread_input%propint(1))
       call mat_free(orbspread_input%propint(6))
       call mat_free(orbspread_input%propint(7))
       call mat_free(orbspread_input%propint(9))
       
       ! propint(5)   -> SECX+SECY+SECZ 
       call mat_daxpy(1E0_realk,orbspread_input%propint(8),orbspread_input%propint(5))
       call mat_free(orbspread_input%propint(8))
       call mat_daxpy(1E0_realk,orbspread_input%propint(10),orbspread_input%propint(5))
       call mat_free(orbspread_input%propint(10))
    ELSE !memory optimized version
       I = 5
       call mat_init(orbspread_input%propint(I),n,n)
       call II_get_single_carmom(6,6,ls%setting,orbspread_input%propint(I),I,nderiv,0E0_realk,0E0_realk,0E0_realk)
       I = 8
       call mat_init(orbspread_input%propint(I),n,n)
       call II_get_single_carmom(6,6,ls%setting,orbspread_input%propint(I),I,nderiv,0E0_realk,0E0_realk,0E0_realk)
       ! propint(5)   -> SECX+SECY+SECZ 
       call mat_daxpy(1E0_realk,orbspread_input%propint(8),orbspread_input%propint(5))
       call mat_free(orbspread_input%propint(8))
       I = 10
       call mat_init(orbspread_input%propint(I),n,n)
       call II_get_single_carmom(6,6,ls%setting,orbspread_input%propint(I),I,nderiv,0E0_realk,0E0_realk,0E0_realk)
       ! propint(5)   -> SECX+SECY+SECZ 
       call mat_daxpy(1E0_realk,orbspread_input%propint(10),orbspread_input%propint(5))
       call mat_free(orbspread_input%propint(10))                
       do I = 2,4
          call mat_init(orbspread_input%propint(I),n,n)
          call II_get_single_carmom(6,6,ls%setting,orbspread_input%propint(I),I,nderiv,0E0_realk,0E0_realk,0E0_realk)
       enddo
    ENDIF
    
  end subroutine orbspread_propint



  subroutine orbspread_init(orbspread_input,m,norb,DoNotAllocateTmpM)
    implicit none
    type(orbspread_data), intent(inout) :: orbspread_input
    integer             , intent(in) :: norb, m
    logical,optional :: DoNotAllocateTmpM
    TYPE(lsitem) :: ls
    integer   :: i
    logical :: AllocTMP

    ! init R Q and tmpM matrices
    call mat_init(orbspread_input%Q,norb,norb)
    do i=1,3
       call mat_init(orbspread_input%R(i),norb,norb)
    enddo

    AllocTMP = .TRUE.
    IF(present(DoNotAllocateTmpM))THEN
       AllocTMP = .NOT.DoNotAllocateTmpM
    ENDIF
    
    IF(AllocTMP)THEN
!       do i=1,4
!          call mat_init(orbspread_input%tmpM(i),norb,norb)
!       enddo
    ENDIF

    ! set norb
    orbspread_input%norb =  norb

    ! set m power
    orbspread_input%m = m

    ! allocate spread2

    call mem_alloc(orbspread_input%spread2,norb)

  end subroutine orbspread_init


  subroutine orbspread_update(orbspread_input,CMO)
    implicit none
    type(orbspread_data), intent(inout) :: orbspread_input
    type(Matrix), intent(in) :: CMO

    type(Matrix) :: tmp
    integer :: nbas,norb,i
    real(realk), pointer :: tmpv(:)

    nbas=CMO%nrow
    norb=CMO%ncol

    call mat_init(tmp,nbas,norb)

    !   R(1:3)
    do i=1,3
       call mat_mul(orbspread_input%propint(i+1),CMO,'n','n',1E0_realk,0E0_realk,tmp)
       call mat_mul(CMO,tmp,'t','n',1E0_realk,0E0_realk,orbspread_input%R(i))
    enddo

    !   Q
    call mat_mul(orbspread_input%propint(5),CMO,'n','n',1E0_realk,0E0_realk,tmp)
    call mat_mul(CMO,tmp,'t','n',1E0_realk,0E0_realk,orbspread_input%Q)

    call mat_free(tmp)

    !   spread2
    call mem_alloc(tmpv,norb)

    call mat_extract_diagonal(orbspread_input%spread2,orbspread_input%Q)

    do i=1,3
       call mat_extract_diagonal(tmpv,orbspread_input%R(i))
       tmpv=tmpv**2
       call daxpy(norb,-1E0_realk,tmpv,1,orbspread_input%spread2,1)
    enddo

    !write(6,*) sqrt(orbspread_input%spread2)
    call mem_dealloc(tmpv)

  end subroutine orbspread_update

  subroutine orbspread_value(oVal, inp)
    real(realk), intent(out) :: oVal
    type(orbspread_data), intent(in)  :: inp

    integer  :: i

    oVal = 0E0_realk
    do i=1,inp%Q%ncol
       oVal = oVal + inp%spread2(i)**inp%m
    enddo


  end subroutine orbspread_value

  subroutine orbspread_update_trustradius(arh,r)
    type(solverItem), intent(inout)   :: arh
    real(realk)     , intent(in)      :: r

    real(realk)                       :: nom, denom

    !nom =   (oVal - arh%old_energy)
    !denom=  arh%denom

    !r = nom/denom 

    if (r.lt. 0E0_realk) then 
       arh%set_max_element = arh%set_max_element/2.0E0_realk
       arh%set_max_step = arh%set_max_step/2.0E0_realk
       arh%step_accepted = .false.
       write (6,*) 'Reject and contract *0.5'
       return
    endif

    if(r.gt. 0.75E0_realk) then
       arh%set_max_element = arh%set_max_element*1.2E0_realk
       arh%set_max_step = arh%set_max_step*1.2E0_realk
       write(6,*) 'Expand *1.2', r
    else if(r.gt. 0.25) then
       write(6,*) 'Keep', r
       continue
    else 
       arh%set_max_element = arh%set_max_element*0.7E0_realk
       arh%set_max_step = arh%set_max_step*0.7E0_realk
       write(6,*) 'Contract *0.7', r
    endif

    arh%step_accepted = .true.

  end subroutine orbspread_update_trustradius




  subroutine orbspread_gradx(G,norb,inp)
    implicit none
    Type(Matrix) , target, intent(inout) :: G
    integer                              :: norb
    type(orbspread_data), intent(inout)  :: inp

    real(realk)  :: diagR(norb,3), tmp(norb)
    integer      :: x, m,i
    type(matrix) :: tmpM1

    m=inp%m

    do x=1, 3
       call mat_extract_diagonal(diagR(:,x),inp%R(x))
    enddo

    tmp =  (inp%spread2**(m-1))
    call mat_zero(G)
    call mat_dmul(tmp,inp%Q,'n',-2E0_realk*m,0E0_realk,G)

    do x=1,3
       do i=1,norb
          tmp(i) = diagR(i,x)*(inp%spread2(i)**(m-1))
       enddo
       call mat_dmul(tmp,inp%R(x),'n',4E0_realk*m,1E0_realk,G)
    enddo

    call mat_init(tmpM1,G%ncol,G%nrow)
    call mat_trans(G,tmpM1)
    call mat_daxpy(-1E0_realk,tmpM1,G)
    call mat_free(tmpM1)

    !call mat_scal(0.5E0_realk,G)

    inp%G => G

  end subroutine  orbspread_gradx
  
!> \brief Compute diagonal hessian elements for variance localization
!> \author Ida-Marie Hoeyvik
  subroutine orbspread_precond_matrix2(inp,P,norb)
    implicit none
    type(orbspread_data), intent(in) :: inp
    type(matrix) :: P
    integer :: m,norb
    real(realk) :: diagQ(norb),diagR(norb,3)
    real(realk) :: spm1(norb),spm2(norb),tmpV(norb)
    integer :: i,x
    type(matrix) :: tmp,tmp1
    real(realk),pointer :: tmpvec(:)

    call mat_zero(P)
    m=inp%m

    call mat_extract_diagonal(diagQ,inp%Q)
    do x=1,3
       call mat_extract_diagonal(diagR(:,x),inp%R(x))
    end do

    do i=1,norb
       spm1(i) = inp%spread2(i)**(m-1)
    enddo
    if (m>1)then
       do i=1,norb
          spm2(i)=inp%spread2(i)**(m-2)
       enddo
    endif
    call mat_dger(dble(2*m),diagQ,spm1,P)

    do x=1,3
       do i=1,norb
          tmpV(i) = spm1(i)*diagR(i,x)
       enddo
       call mat_dger(-dble(4*m),tmpV,diagR(:,x),P)
    end do

    call mat_init(tmp,P%nrow,P%ncol)
    call mat_zero(tmp)
    if (m>1) then
       call mat_hmul(1d0,inp%Q,inp%Q,0d0,tmp)
       call mat_dmul(spm2,tmp,'n',dble(4*m*(m-1)),1d0,P)
    end if
    do x=1,3
       call mat_hmul(1d0,inp%R(x),inp%R(x),0d0,tmp)
       call mat_dmul(spm1,tmp,'n',-dble(8*m),1d0,P)
       if (m>1) then
          call mat_hmul(1d0,inp%Q,inp%R(x),0d0,tmp)
          call mat_dmul(spm2*diagR(:,x),tmp,'n',-dble(16*m*(m-1)),1d0,P)
       end if
    end do

    call mat_zero(tmp)
    if (m>1) then
       do x=1,3
          call mat_dmul(diagR(:,x),inp%R(x),'n',1d0,1d0,tmp)
       end do
       call mat_init(tmp1,P%nrow,P%ncol)
       call mat_dmul(spm2,tmp,'n',dble(16*m*(m-1)),0d0,tmp1)
       call mat_hmul(1d0,tmp,tmp1,1d0,P)
       call mat_free(tmp1)
    end if

    call mem_alloc(tmpvec,norb)
    tmpvec=1.0d0
    call mat_dger(-dble(2*m),diagQ*spm1,tmpvec,P)

    do x=1,3
       do i=1,norb
          tmpV(i) = spm1(i)*diagR(i,x)*diagR(i,x)
       enddo
       call mat_dger(dble(4*m),tmpV,tmpvec,P)
    end do
    call mem_dealloc(tmpvec)

    call mat_trans(P,tmp)
    call mat_daxpy(1d0,tmp,P)
    call mat_free(tmp)
    call mat_scal_dia(0d0,P)

  end subroutine orbspread_precond_matrix2

end module orbspread_utilMod
