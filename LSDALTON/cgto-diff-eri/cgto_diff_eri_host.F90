MODULE cgto_diff_eri_host_interface
  use precision
  use molecule_typeType, only: MOLECULEINFO
  use AO_TypeType, only: AOITEM,maxAOangmom
  use AO_Type, only: free_aoitem
  use integral_type, only: INTEGRALINPUT
  use TYPEDEFTYPE, only: LSSETTING
  use BUILDAOBATCH, only: BUILD_AO
  use memory_handling
#ifdef BUILD_CGTODIFF
  use interface_basis !cgto_diff_eri code
#endif
  use Matrix_module, only: Matrix
  use Matrix_operations
contains
  subroutine cgto_diff_eri_DGD_Econt(D,Econt,ls_setting, luprint)
    implicit none
    integer,intent(in) :: luprint
    type(matrix),intent(in)    :: D
    type(LSSETTING), intent(inout) :: ls_setting
    real(realk), intent(inout)  :: Econt
    !
#ifdef BUILD_CGTODIFF
    type(AOITEM) :: ao_items  
    integer :: istart,istartCC,I,J,nelms,nEXPSIZE,nCCSIZE,nbast,icc
    logical :: nofamily
    integer,pointer :: ICHARGE(:),IBASIS(:),CCINDEX(:),ANGMOM(:),AMOM(:)
    integer,pointer :: nPrim(:),nCont(:),start_exponents(:),start_CC(:)
    real(realk),pointer :: CENTER(:,:),exponents(:),ContC(:),average(:)
    real(realk),pointer :: Dfull(:,:)
    !==========================================================
    ! Step 1: init the basis 
    !==========================================================
    nbast = D%nrow
    nofamily = ls_setting%SCHEME%nofamily
    ls_setting%SCHEME%nofamily = .TRUE.
    call BUILD_AO(luprint,ls_setting%SCHEME,ls_setting%SCHEME%AOPRINT,&
         & ls_setting%molecule(1)%p,ls_setting%Basis(1)%p%REGULAR, ao_items,&
         & .false.,.false.,.false.)
    
    call mem_alloc(CENTER,3,AO_items%nbatches)
    call mem_alloc(ANGMOM,AO_items%nbatches)
    call mem_alloc(CCINDEX,AO_items%nbatches)
    call mem_alloc(ICHARGE,AO_items%nbatches)
    call mem_alloc(IBASIS,AO_items%nbatches)
    IF(AO_items%nCC.NE.AO_items%nEXP)THEN
       print*,'AO_items%nCC',AO_items%nCC
       print*,'AO_items%nEXP',AO_items%nEXP
       CALL LSQUIT('AO_items%nCC.NE.AO_items%nEXP',-1)
    ENDIF
    do I=1,AO_items%nbatches
       CENTER(1,I) = AO_items%BATCH(I)%CENTER(1)
       CENTER(2,I) = AO_items%BATCH(I)%CENTER(2)
       CENTER(3,I) = AO_items%BATCH(I)%CENTER(3)
       ANGMOM(I) = AO_items%BATCH(I)%ANGMOM(1)
       IF(AO_items%BATCH(I)%nANGMOM.GT.1)call lsquit('use .NOFAMILY',-1)
       do J=1,AO_items%BATCH(I)%pExponents%nrow
!          exponents(istart+J) = AO_items%BATCH(I)%pExponents%elms(J)
          IF(ABS((AO_items%BATCH(I)%pExponents%elms(J)-AO_items%Exponents(AO_items%BATCH(I)%CCindex(1))%elms(J))).GT.1.0E-10_realk)THEN
             CALL LSQUIT('pExponents%elms(J)',-1)
          ENDIF
       enddo       
       CCINDEX(I) = AO_items%BATCH(I)%CCindex(1)

       nelms = AO_items%BATCH(I)%pCC(1)%p%nrow*AO_items%BATCH(I)%pCC(1)%p%ncol
       do J=1,nelms
!          CC(istartCC+J) = AO_items%BATCH(I)%pCC%elms(J)
          IF(ABS(AO_items%BATCH(I)%pCC(1)%p%elms(J)-AO_items%CC(AO_items%BATCH(I)%CCindex(1))%elms(J)).GT.1.0E-10_realk)THEN
             CALL LSQUIT('pCC%elms(J)',-1)
          ENDIF
       enddo
       ICHARGE(I) = HUGE(nelms) !ajt FIXME wrong
       IBASIS(I)=I-1
    enddo
    
    call mem_alloc(nprim,AO_items%nCC)
    call mem_alloc(nCont,AO_items%nCC)
    call mem_alloc(start_exponents,AO_items%nCC)
    call mem_alloc(start_CC,AO_items%nCC)
    istart = 0
    istartCC = 0
    do I=1,AO_items%nCC
       start_exponents(I) = istart !start from zero
       istart = istart + AO_items%CC(I)%nrow
       nprim(I) = AO_items%CC(I)%nrow
       nCont(I) = AO_items%CC(I)%ncol
       nelms = AO_items%CC(I)%nrow*AO_items%CC(I)%ncol
       start_CC(I) = istartCC !start from zero
       istartCC = istartCC + nelms
    enddo

    call mem_alloc(exponents,istart)
    call mem_alloc(ContC,istartCC)
    nEXPSIZE = istart
    nCCSIZE = istartCC

    istart = 0
    istartCC = 0
    do I=1,AO_items%nCC
       do J=1,AO_items%Exponents(I)%nrow
          exponents(istart+J) = AO_items%Exponents(I)%elms(J)
       enddo
       istart = istart + AO_items%CC(I)%nrow
       
       nelms = AO_items%CC(I)%nrow*AO_items%CC(I)%ncol
       do J=1,nelms
          ContC(istartCC+J) = AO_items%CC(I)%elms(J)
       enddo
       istartCC = istartCC + nelms
    enddo

    !build amom    
    call mem_alloc(amom,AO_items%nCC)
    do ICC=1,AO_items%nCC
       amom(icc) = -1
    enddo
    do I=1,AO_items%nbatches
       icc = CCINDEX(I)
       amom(icc) = ANGMOM(I)
    enddo
    print*,'interface_basis_init_general'
    !builds the cgto basis in (module interface_basis)
    !external/openrsp/src/interfaces/interface_basis.F90 
    call interface_basis_init_general(nbast,AO_items%nbatches,&
         & nEXPSIZE,nCCSIZE,AO_items%nCC,ANGMOM,CCINDEX,ICHARGE,IBASIS,nPrim,&
         & nCont,start_exponents,start_CC,CENTER,exponents,ContC,&
         & maxAOangmom,amom)
    call mem_alloc(average,1)

    call mem_alloc(Dfull,nbast,nbast)
    call mat_to_full(D, 1.0E0_realk, Dfull)
    average(1) = 0.0E0_realk
    call interface_basis_init_argument_general(Dfull,Average,1,0,nbast)

    call interface_basis_main_general()
    print *, '2-el energy from cgto-diff-eri', average(1)
    Econt = average(1)
    call interface_basis_finalize()

    !Free stuff
    call mem_dealloc(Dfull)
    call FREE_AOITEM(luprint,AO_items)
    call mem_dealloc(amom)
    call mem_dealloc(CENTER)
    call mem_dealloc(ANGMOM)
    call mem_dealloc(CCINDEX)
    call mem_dealloc(ICHARGE)
    call mem_dealloc(IBASIS)
    call mem_dealloc(nprim)
    call mem_dealloc(nCont)
    call mem_dealloc(start_exponents)
    call mem_dealloc(start_CC)
    call mem_dealloc(exponents)
    call mem_dealloc(ContC)
    call mem_dealloc(average)

    !restore settings
    ls_setting%SCHEME%nofamily = nofamily
#endif
  end subroutine cgto_diff_eri_DGD_Econt

  subroutine cgto_diff_eri_xfac_general(ExchangeFactor)
    implicit none
    real(realk) :: ExchangeFactor
#ifdef BUILD_CGTODIFF
    call interface_basis_xfac_general(ExchangeFactor)
#endif
  end subroutine cgto_diff_eri_xfac_general

END MODULE CGTO_DIFF_ERI_HOST_INTERFACE
