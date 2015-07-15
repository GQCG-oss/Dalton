!> @file
!> Simple integrals related
!> \author Marcin Ziolkowski, (Pablo Baudin: added subroutines to get 
!          (non-T1-transformed MO int. used for MO-based CCSD)
module ccintegrals

  use memory_handling
  use precision
  use lstiming!, only: lstimer
  use dec_typedef_module
  use typedeftype!, only: Lsitem
  use LSTIMING!,only:lstimer
  use matrix_module!, only:matrix
  use matrix_operations
  use integralinterfaceDEC
  use screen_mod
  use LSparameters
  use integralinterfaceMOD
  use II_XC_interfaceModule
  use IchorErimoduleHost

  ! MO-CCSD module:
  use tensor_interface_module
  use buildaobatch
  use memory_handling
  use daltoninfo
  use lspdm_tensor_operations_module
#ifdef VAR_MPI
  use lsmpi_type
  use infpar_module
#endif

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use array2_simple_operations
  use array4_simple_operations
  ! MO-CCSD module:
  use dec_fragment_utils 
  use cc_tools_module
#ifdef VAR_MPI
  use decmpi_module
#endif

  interface getL
     module procedure getL_simple
     module procedure getL_simple_from_gmo
     module procedure getL_diff
  end interface

  private :: get_mem_t1_free_gmo, get_mem_MO_CCSD_residual, &
       & init_gmo_arrays, gao_to_gmo, get_MO_batches_info, pack_and_add_gmo

contains

  !> \brief Get full two-electron integrals in AO basis using DALTONs integral program
  subroutine get_full_eri(mylsitem,nbasis,g_ao)

    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer, intent(in) :: nbasis
    type(array4),intent(inout) :: g_ao
    type(matrix) :: g_matrix
    integer, dimension(4) :: ao_dims
    character :: intspec(5)
    logical :: SAMEMOL
    integer :: iprint
    intspec(1) = 'R'
    intspec(2) = 'R'
    intspec(3) = 'R'
    intspec(4) = 'R'
    intspec(5) = 'C'
    ao_dims=[nbasis,nbasis,nbasis,nbasis]

    g_ao = array4_init_standard(ao_dims)
    ! KK Quick fix: Filename associated with g_ao
    g_ao%filename = 'gao'
    IF(DECinfo%useIchor)THEN
       SameMOL = .TRUE.
       iprint = 0
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,&
            & INTSPEC,SameMOL)
       call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,&
            & nbasis,nbasis,nbasis,nbasis,g_ao%val,INTSPEC,.TRUE.,&
            & 1,1,1,1,1,1,1,1,.FALSE.,nbasis,nbasis,nbasis,nbasis,.FALSE.,&
            & DECinfo%IntegralThreshold)
       call FREE_SCREEN_ICHORERI()
    ELSE
       call ii_get_4center_eri(DECinfo%output,DECinfo%output,&
            & mylsitem%setting,g_ao%val,nbasis,nbasis,nbasis,nbasis,intspec)
    ENDIF
    ! write ao integrals to disk
    call array4_write(g_ao)

  end subroutine get_full_eri


  !> \brief Get L as L[p,q,r,s,]=2g[p,q,r,s]-g[p,s,r,q]
  function getL_simple(gao,left1,left2,right1,right2) result(l)

    implicit none
    type(array2), intent(in) :: left1,left2,right1,right2
    type(array4), intent(inout) :: gao
    type(array4) :: tmp
    type(array4) :: l

    l = get_gmo_simple(gao,left1,left2,right1,right2)
    tmp = array4_duplicate(l)
    call array4_reorder(tmp,[1,4,3,2])
    call array4_scale(l,2.0E0_realk)
    call array4_add_to(l,-1.0E0_realk,tmp)
    call array4_free(tmp)

    return
  end function getL_simple

  !> \brief Get L from gmo
  function getL_simple_from_gmo(gmo) result(l)

    implicit none
    type(array4) :: gmo ! Not changed
    type(array4) :: l

    if(DECinfo%array4OnFile) then ! array values stored on file
       l = getL_simple_from_gmo_file(gmo)
    else ! array values stored in memory
       l = getL_simple_from_gmo_memory(gmo)
    end if

    return
  end function getL_simple_from_gmo


  !> \brief Get L from gmo, storing values in memory.
  function getL_simple_from_gmo_memory(gmo) result(l)

    implicit none
    type(array4) :: gmo ! Not changed
    type(array4) :: l

    l = array4_duplicate(gmo)
    call array4_scale(l,2.0E0_realk)
    call array4_reorder(gmo,[1,4,3,2])
    call array4_add_to(l,-1.0E0_realk,gmo)
    call array4_reorder(gmo,[1,4,3,2])

    return
  end function getL_simple_from_gmo_memory



  !> \brief Get L from gmo, storing values on file.
  !> L_{bjai} = 2*g_{bjai} - g_{ajbi}
  !> \author Kasper Kristensen
  !> \date October 2010
  function getL_simple_from_gmo_file(gmo) result(l)

    implicit none
    type(array4) :: gmo ! Not changed
    type(array4) :: l
    integer :: i,j,a,b,c,d,nocc,nvirt
    real(realk),pointer :: gtmp(:,:,:), Ltmp(:,:,:)

    ! gmo stored as (virt,occ,virt,occ)
    nvirt = gmo%dims(1)
    nocc  = gmo%dims(2)


    ! Sanity checks
    ! *************

    ! Here we assume that the values in gmo are stored using storing type 2.
    if(gmo%storing_type /= 2) then
       call lsquit('getL_simple_from_gmo_file: &
            & Only implemented when for storing type 2!', DECinfo%output)
    end if

    ! gmo stored as (virt,occ,virt,occ)
    if( (gmo%dims(3) /= nvirt) .or.  (gmo%dims(4) /= nocc) ) then
       call lsquit('getL_simple_from_gmo_file: &
            & gmo dimensions must be (virt,occ,virt,occ)',DECinfo%output)
    end if


    ! Initialize stuff
    ! ****************

    ! Temporary arrays for reading and writing
    call mem_alloc(gtmp,nvirt,nocc,nvirt)
    call mem_alloc(Ltmp,nvirt,nocc,nvirt)

    ! Initialize L using storing type 2
    L = array4_init([nvirt,nocc,nvirt,nocc],2,.false.)

    ! Open files
    call array4_open_file(gmo)
    call array4_open_file(L)


    ! Calculate: L_{bjai} = 2*g_{bjai} - g_{ajbi}
    ! *******************************************

    i_loop: do i=1,nocc

       c_loop: do c=1,nvirt
          ! Read in gmo(:,:,c,i)
          call array4_read_file(gmo,c,i,gtmp(:,:,c),nvirt,nocc)
       end do c_loop

       ! Now gtmp contains gmo(:,:,:,i) and L(:,:,:,i) can be constructed

       a_loop: do a=1,nvirt
          j_loop: do j=1,nocc
             b_loop: do b=1,nvirt
                ! Calculate L_{b,j,a,i}
                Ltmp(b,j,a) = 2E0_realk*gtmp(b,j,a) - gtmp(a,j,b)
             end do b_loop
          end do j_loop
       end do a_loop

       ! At this point L(:,:,:,i) has been calculated and is stored in Ltmp
       ! Write Ltmp to file using storing type 2
       ! I.e. for each i we write Ltmp(:,:,d) to L(:,:,d,i)
       d_loop: do d=1,nvirt
          call array4_write_file(L,d,i,Ltmp(:,:,d),nvirt,nocc)
       end do d_loop

    end do i_loop



    ! Free stuff
    call array4_close_file(gmo,'keep')
    call array4_close_file(L,'keep')
    call mem_dealloc(gtmp)
    call mem_dealloc(Ltmp)

  end function getL_simple_from_gmo_file



  !> \brief Get L array
  function getL_diff(gao,A_left1,A_left2,A_right1,A_right2, &
       B_left1,B_left2,B_right1,B_right2) result(l)

    implicit none
    type(array4), intent(inout) :: gao
    type(array2), intent(in) :: A_left1,A_left2,A_right1,A_right2, &
         B_left1,B_left2,B_right1,B_right2
    type(array4) :: tmp
    type(array4) :: l
    integer :: p,q,r,s

    l = get_gmo_simple(gao,A_left1,A_left2,A_right1,A_right2)
    tmp = get_gmo_simple(gao,B_left1,B_left2,B_right1,B_right2)

    !#ifdef EXTRA_SIMPLE_INTEGRALS
    do p = 1,l%dims(1)
       do q = 1,l%dims(2)
          do r = 1,l%dims(3)
             do s = 1,l%dims(4)
                l%val(p,q,r,s) = 2.0E0_realk*l%val(p,q,r,s) - tmp%val(p,s,r,q)
             end do
          end do
       end do
    end do
    !#else
    !    call array4_scale(l,2.0E0_realk)
    !    call array4_reorder(l,[2,4,1,3]) ! l[pq,rs] -> l[qs,pr]
    !    call array4_reorder(tmp,[2,4,1,3]) ! tmp[ps,rq] -> tmp[sq,pr]
    !    call array4_add_to(l,-1.0E0_realk,tmp)
    !    call array4_reorder(l,[3,1,4,2])
    !#endif

    call array4_free(tmp)

    return
  end function getL_diff

  !> \brief Simple routine to transform integrals to mo (general MO integrals)
  function get_gmo_simple(gao,left1,left2,right1,right2) result(gmo)

    implicit none
    type(array4), intent(inout) :: gao
    type(array4) :: gmo
    type(array4) :: tmp1,tmp2,tmp3
    type(array2), intent(in) :: left1,left2,right1,right2
    real(realk), pointer :: AA(:,:), BB(:,:), CC(:,:), DD(:,:)
    real(realk) :: starttime,endtime
    integer, dimension(4) :: dims,rotate
    integer :: nbasis
    integer :: a,b,c,d,p,q,r,s

    nbasis = left1%dims(1)

    dims(1) = left1%dims(2)
    dims(2) = left2%dims(2)
    dims(3) = right1%dims(2)
    dims(4) = right2%dims(2)

#ifdef EXTRA_SIMPLE_INTEGRALS
    call array4_read(gao)

    AA => left1%val
    BB => left2%val
    CC => right1%val
    DD => right2%val

    gmo = array4_init(dims)

    tmp1 = array4_init([left1%dims(2),nbasis,nbasis,nbasis])
    tmp2 = array4_init([left1%dims(2),left2%dims(2),nbasis,nbasis])
    tmp3 = array4_init([left1%dims(2),left2%dims(2),right1%dims(2),nbasis])

    do d=1,nbasis
       do c=1,nbasis
          do b=1,nbasis
             do a=1,left1%dims(2)

                do p=1,nbasis
                   tmp1%val(a,b,c,d) = tmp1%val(a,b,c,d) + &
                        + gao%val(p,b,c,d)*AA(p,a)
                end do

             end do
          end do
       end do
    end do

    print *,'trans 1 done'

    do d=1,nbasis
       do c=1,nbasis
          do b=1,left2%dims(2)
             do a=1,left1%dims(2)

                do p=1,nbasis
                   tmp2%val(a,b,c,d) = tmp2%val(a,b,c,d) + tmp1%val(a,p,c,d)*BB(p,b)
                end do

             end do
          end do
       end do
    end do

    print *,'trans 2 done'

    do d=1,nbasis
       do c=1,right1%dims(2)
          do b=1,left2%dims(2)
             do a=1,left1%dims(2)

                do p=1,nbasis
                   tmp3%val(a,b,c,d) = tmp3%val(a,b,c,d) + tmp2%val(a,b,p,d)*CC(p,c)
                end do

             end do
          end do
       end do
    end do

    print *,'trans 3 done'

    do d=1,right2%dims(2)
       do c=1,right1%dims(2)
          do b=1,left2%dims(2)
             do a=1,left1%dims(2)

                do p=1,nbasis
                   gmo%val(a,b,c,d) = gmo%val(a,b,c,d) + tmp3%val(a,b,c,p)*DD(p,d)
                end do

             end do
          end do
       end do
    end do

    print *,'trans 4 done'

    call array4_free(tmp1)
    call array4_free(tmp2)
    call array4_free(tmp3)

    AA => null()
    BB => null()
    CC => null()
    DD => null()

    call array4_dealloc(gao)

#else

    call cpu_time(starttime)

    rotate=[2,3,4,1]

    ! read ao integrals
    call array4_read(gao)
    tmp1 = array4_init([left1%dims(2),nbasis,nbasis,nbasis])
    call array4_contract1(gao,left1,tmp1,.true.)

    ! deallocate ao integrals
    call array4_dealloc(gao)

    tmp2 = array4_init([left2%dims(2),nbasis,nbasis,left1%dims(2)])
    call array4_reorder(tmp1,rotate)
    call array4_contract1(tmp1,left2,tmp2,.true.)
    call array4_free(tmp1)

    tmp3 = array4_init([right1%dims(2),nbasis,left1%dims(2),left2%dims(2)])
    call array4_reorder(tmp2,rotate)
    call array4_contract1(tmp2,right1,tmp3,.true.)
    call array4_free(tmp2)

    gmo = array4_init([right2%dims(2),left1%dims(2),left2%dims(2),right1%dims(2)])
    call array4_reorder(tmp3,rotate)
    call array4_contract1(tmp3,right2,gmo,.true.)
    call array4_reorder(gmo,rotate)
    call array4_free(tmp3)

    call cpu_time(endtime)

#endif

    return
  end function get_gmo_simple

  !> \brief Transform two indices (to construct Fock matrix)
  function get_oopq(gao,left1,left2) result(gmo)

    implicit none
    type(array4), intent(inout) :: gao
    type(array4) :: gmo
    type(array4) :: tmp
    type(array2), intent(inout) :: left1,left2
    integer, dimension(4) :: dims
    integer :: nbasis

    nbasis = left1%dims(1)

    dims(1) = left1%dims(2)
    dims(2) = left2%dims(2)
    dims(3) = nbasis
    dims(4) = nbasis

    ! read ao integrals
    call array4_read(gao)
    tmp = array4_init([left2%dims(2),nbasis,nbasis,nbasis])
    call array4_contract1(gao,left2,tmp,.true.)

    ! deallocate ao integrals
    call array4_dealloc(gao)

    gmo = array4_init([left1%dims(2),left2%dims(2),nbasis,nbasis])
    call array4_reorder(tmp,[2,1,3,4])
    call array4_contract1(tmp,left1,gmo,.true.)
    call array4_free(tmp)

    return
  end function get_oopq

  !> \brief Transform two indices (to construct Fock matrix)
  function get_exchange_as_oopq(gao,left1,left2) result(gmo)

    implicit none
    type(array4), intent(inout) :: gao
    type(array4) :: gmo
    type(array4) :: tmp
    type(array2), intent(inout) :: left1,left2
    integer, dimension(4) :: dims
    integer :: nbasis

    nbasis = left1%dims(1)

    dims(1) = left1%dims(2)
    dims(2) = left2%dims(2)
    dims(3) = nbasis
    dims(4) = nbasis

    ! read ao integrals
    call array4_read(gao)
    call array4_reorder(gao,[3,2,1,4])
    tmp = array4_init([left2%dims(2),nbasis,nbasis,nbasis])
    call array4_contract1(gao,left2,tmp,.true.)

    ! deallocate ao integrals
    call array4_reorder(gao,[3,2,1,4])
    call array4_dealloc(gao)

    gmo = array4_init([left1%dims(2),left1%dims(2),nbasis,nbasis])
    call array4_reorder(tmp,[2,1,3,4])
    call array4_contract1(tmp,left1,gmo,.true.)
    call array4_free(tmp)

    return
  end function get_exchange_as_oopq


  !> \brief Wrapper for dec_fock_transformation using fortran arrays, only to be used
  !> for testing purposes!
  !> \author Kasper Kristensen
  !> \date June 2015
  subroutine dec_fock_transformation_fortran_array(nbasis,FockU,U,MyLsitem,symmetric,incl_h)

    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Fock transformation on U matrix
    real(realk), intent(inout) :: FockU(nbasis,nbasis)
    !> Matrix to carry out Fock transformation on
    real(realk), intent(inout) :: U(nbasis,nbasis)
    !> LS DALTON info
    type(lsitem), intent(inout) :: MyLsItem
    !> Is U symmetric (true) or not (false)?
    logical, intent(in) :: symmetric
    !> Include one-electron contribution to Fock matrix (default: not include)
    logical,intent(in),optional :: incl_h
    type(matrix) :: FockU_mat, U_mat

    call mat_init(FockU_mat,nbasis,nbasis)
    call mat_init(U_mat,nbasis,nbasis)
    call mat_set_from_full(FockU,1E0_realk, FockU_mat)
    call mat_set_from_full(U,1E0_realk, U_mat)
    if(present(incl_h)) then
       call dec_fock_transformation(FockU_mat,U_mat,MyLsitem,symmetric,incl_h=incl_h)
    else
       call dec_fock_transformation(FockU_mat,U_mat,MyLsitem,symmetric)
    end if
    call mat_to_full(FockU_mat, 1.0_realk, FockU)
    call mat_to_full(U_mat, 1.0_realk, U)

    call mat_free(FockU_mat)
    call mat_free(U_mat)

  end subroutine dec_fock_transformation_fortran_array


  !> \brief Carry out 2-electron Fock transformation on matrix U, i.e.
  !> FockU = 2J(U) - K(U)
  !> where J(U) and K(U) are the coulomb and exchange transformations, respectively.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_fock_transformation(FockU,U,MyLsitem,symmetric,incl_h)

    implicit none
    !> Fock transformation on U matrix
    type(matrix), intent(inout) :: FockU
    !> Matrix to carry out Fock transformation on
    type(matrix), intent(in) :: U
    !> LS DALTON info
    type(lsitem), intent(inout) :: MyLsItem
    !> Is U symmetric (true) or not (false)?
    logical, intent(in) :: symmetric
    !> Include one-electron contribution to Fock matrix (default: not include)
    logical,intent(in),optional :: incl_h
    real(realk) :: Edft(1),DFTELS
    logical :: doMPI 
    integer :: igrid
    type(matrix) :: h
    ! Sanity check
    if(U%nrow /= U%ncol) then
       call lsquit('dec_fock_transformation:&
            & Matrix U must be quadratic!',-1)
    end if

    ! Carry out Fock transformation on U
    call II_get_Fock_mat(DECinfo%output, DECinfo%output, &
         & MyLsitem%setting,U,symmetric,FockU,1,.FALSE.)
    IF(DECinfo%DFTreference)THEN
       !rebuild the grid 
       MyLsitem%setting%scheme%DFT%griddone = 0
       igrid = MyLsitem%setting%scheme%DFT%igrid
       MyLsitem%setting%scheme%DFT%GridObject(igrid)%GRIDDONE = 0

       !not the full number of electrons we deactivate the testing
       DFTELS = MyLsItem%setting%scheme%DFT%DFTELS
       MyLsItem%setting%scheme%DFT%DFTELS = 100.0E0_realk

       !Deactivate MPI - done at the DEC level - maybe this could
       !be done in the local group - then the MyLsItem%setting%node
       !and MyLsItem%setting%comm needs to be set correctly
       !       doMPI = MyLsItem%setting%scheme%doMPI
       !       MyLsItem%setting%scheme%doMPI = .FALSE.
       print*,'MyLsItem%setting%node',MyLsItem%setting%node
       print*,'MyLsItem%setting%numnodes',MyLsItem%setting%numnodes
       call II_get_xc_fock_mat(DECinfo%output,DECinfo%output,&
            & MyLsItem%setting,U%nrow,U,FockU,Edft,1)       
       !       MyLsItem%setting%scheme%doMPI = doMPI 
       !       IF(DECinfo%FrozenCore)
       MyLsItem%setting%scheme%DFT%DFTELS = DFTELS
    ENDIF


    ! Add one-electron contribution
    if(present(incl_h)) then
       if(incl_h) then
          call mat_init(h,U%nrow,U%ncol)
          call mat_zero(h)
          call ii_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,h)
          call mat_daxpy(1E0_realk,h,FockU)
          call mat_free(h)
       end if
    end if

  end subroutine dec_fock_transformation

  !> \brief Calculate full AO integrals. Only for testing.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine get_full_AO_integrals(nbasis,ncabs,gao,MyLsitem,intSpec)

    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis,ncabs
    !> AO integrals
    real(realk),pointer :: gao(:,:,:,:)
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    character(5),intent(IN) :: intSpec
    TYPE(DECscreenITEM)   :: DecScreen
    logical :: doscreen,SameMOL,NoSymmetry,FULLBATCH,MoTrans
    integer :: ndim(4),i,iAO,nAObatches,iprint
    character(1) :: intSpecConvert(5)

    intSpecConvert(1) = intSpec(1:1)
    intSpecConvert(2) = intSpec(2:2)
    intSpecConvert(3) = intSpec(3:3)
    intSpecConvert(4) = intSpec(4:4)
    intSpecConvert(5) = intSpec(5:5)

    DO i=1,4
       IF (intSpecConvert(i).EQ.'R') THEN !Regular AO basis
          ndim(i) = nbasis
       ELSE IF (intSpecConvert(i).EQ.'C') THEN !CABS AO basis
          ndim(i) = ncabs
       ELSE
          CALL LSQUIT('Error in get_full_AO_integrals in specification of AOs',-1)
       ENDIF
    ENDDO
    IF(DECinfo%useIchor)THEN
       !Use Ichor code to calculate Integrals     
       !Calculate Screening integrals 
       SameMOL = .TRUE. !Specifies same molecule on all centers 
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       iAO = 1 !which center? 
       call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
       ! Get AO integrals
       FULLBATCH = .TRUE.   !Full set of AOs?
       NoSymmetry = .FALSE. !Use Permutational symmetry?
       MoTrans = .FALSE.    !Transformation to Molecular Orbital basis?
       iprint = 0           !Printlevel 
       call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,ndim(1),ndim(2),ndim(3),ndim(4),&
            & gao,INTSPEC,FULLBATCH,1,nAObatches,1,nAObatches,1,nAObatches,1,nAObatches,&
            & MoTrans,ndim(1),ndim(2),ndim(3),ndim(4),NoSymmetry,DECinfo%IntegralThreshold)
       !Free screening info
       call FREE_SCREEN_ICHORERI()
    ELSE
       !Use Thermite code to calculate Integrals     
       ! Set integral screening
       doscreen = mylsitem%setting%scheme%cs_screen .OR. mylsitem%setting%scheme%ps_screen
       call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mylsitem%setting,1,1,intspecConvert,DECinfo%IntegralThreshold)
       IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
       IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabRHS
       
       ! Get AO integrals
       call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output,Mylsitem%SETTING,&
            & gao,1,1,ndim(3),ndim(4),ndim(1),ndim(2),ndim(3),ndim(4),.true.,intSpecConvert,DECinfo%IntegralThreshold)
       call free_decscreen(DECSCREEN)
       nullify(mylsitem%setting%LST_GAB_RHS)
       nullify(mylsitem%setting%LST_GAB_LHS)       
    ENDIF

  end subroutine get_full_AO_integrals

  !> \brief Calculate mixed one-electron and Coulomb matrix contributions.
  !> \author S Reine
  !> \date May 2012
  subroutine get_AO_Fock(nbasis,ncabs,F,D,MyLsitem,intSpec)

    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis,ncabs
    !> Mixed one-electron and Coulomb matrix contribution
    TYPE(Matrix),intent(inout) :: F
    !> SCF AO density-matrixc
    TYPE(Matrix),intent(in)    :: D
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Specified the four AOs and the two-electron operator
    character(5),intent(IN) :: intSpec
    !
    TYPE(Matrix) :: K

    call get_AO_hJ(nbasis,ncabs,F,D,MyLsitem,intSpec)
    call mat_init(K,F%nrow,F%ncol)
    call get_AO_K(nbasis,ncabs,K,D,MyLsitem,intSpec)
    call mat_daxpy(1E0_realk,K,F)
    call mat_free(K)
  end subroutine get_AO_Fock


  !> \brief Calculate mixed one-electron and Coulomb matrix contributions.
  !> \author S Reine
  !> \date May 2012
  subroutine get_AO_hJ(nbasis,ncabs,hJ,D,MyLsitem,intSpec)


    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis,ncabs
    !> Mixed one-electron and Coulomb matrix contribution
    TYPE(Matrix),intent(inout) :: hJ
    !> SCF AO density-matrixc
    TYPE(Matrix),intent(in)    :: D
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Specified the four AOs and the two-electron operator
    character(5),intent(IN) :: intSpec
    integer :: dummy
    TYPE(DECscreenITEM)   :: DecScreen
    logical :: doscreen,doMPI
    integer :: ndim(4),i,AO(4),Oper,ndmat
    TYPE(Matrix) :: h,Jarr(1)
    character(1) :: intSpecConvert(5)

    ndmat=1    

    intSpecConvert(1) = intSpec(1:1)
    intSpecConvert(2) = intSpec(2:2)
    intSpecConvert(3) = intSpec(3:3)
    intSpecConvert(4) = intSpec(4:4)
    intSpecConvert(5) = intSpec(5:5)


    DO i=1,4
       IF (intSpecConvert(i).EQ.'R') THEN !Regular AO basis
          ndim(i) = nbasis
          AO(i)   = AORegular
       ELSE IF (intSpecConvert(i).EQ.'C') THEN !CABS AO basis
          ndim(i) = ncabs
          AO(i)   = AOdfCABS
       ELSE
          CALL LSQUIT('Error in get_full_AO_integrals in specification of AOs',-1)
       ENDIF
    ENDDO

    IF ((AO(3).NE.AORegular).OR.(AO(4).NE.AORegular)) &
         &  CALL LSQUIT('Error in get_full_AO_integrals in AO3 or AO4. Both should be AOregular',-1)

    IF (intSpecConvert(5).EQ.'C') THEN
       Oper = coulombOperator
    ELSE
       CALL LSQUIT('Error in get_full_AO_integrals in specification of operator',-1)
    ENDIF

    ! Quick fix because we need to pass an array
    call mat_init(Jarr(1),hJ%nrow,hJ%ncol)
    call mat_zero(Jarr(1))

    call mat_init(h,ndim(1),ndim(2))
    doMPI = MyLsitem%SETTING%scheme%doMPI
    MyLsitem%SETTING%scheme%doMPI = .NOT.(((AO(1).EQ.AOdfCABS).OR.(AO(2).EQ.AOdfCABS)).OR.&
         & ((AO(3).EQ.AOdfCABS).OR.(AO(4).EQ.AOdfCABS)))
    CALL II_get_h1_mixed(DECinfo%output,DECinfo%output,MyLsitem%SETTING,h,AO(1),AO(2))
    CALL II_get_coulomb_mat_mixed(DECinfo%output,DECinfo%output,MyLsitem%SETTING,(/D/),Jarr,ndmat,&
         &                            AO(1),AO(2),AO(3),AO(4),Oper)
    MyLsitem%SETTING%scheme%doMPI = doMPI

    call mat_assign(hJ,Jarr(1))
    call mat_free(Jarr(1))
    CALL mat_daxpy(1E0_realk,h,hJ)
    CALL mat_free(h)


  end subroutine get_AO_hJ

  !> \brief Calculate mixed one-electron and Coulomb matrix contributions.
  !> \author S Reine
  !> \date May 2012
  subroutine get_AO_K(nbasis,ncabs,K,D,MyLsitem,intSpec)

    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis,ncabs
    !> Mixed exchange-matrix contribution
    TYPE(Matrix),intent(inout) :: K
    !> SCF AO density-matrixc
    TYPE(Matrix),intent(in)    :: D
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Specified the four AOs and the two-electron operator
    character(5),intent(IN) :: intSpec
    integer :: dummy
    TYPE(DECscreenITEM)   :: DecScreen
    logical :: doscreen
    integer :: ndim(4),i,AO(4),Oper,ndmat
    TYPE(Matrix) :: h
    character(1) :: intSpecConvert(5)
    logical :: Dsym,doMPI
    type(matrix) :: Karr(1)

    Dsym=.true.
    ndmat=1

    intSpecConvert(1) = intSpec(1:1)
    intSpecConvert(2) = intSpec(2:2)
    intSpecConvert(3) = intSpec(3:3)
    intSpecConvert(4) = intSpec(4:4)
    intSpecConvert(5) = intSpec(5:5)


    DO i=1,4
       IF (intspecConvert(i).EQ.'R') THEN !Regular AO basis
          ndim(i) = nbasis
          AO(i)   = AORegular
       ELSE IF (intspecConvert(i).EQ.'C') THEN !CABS AO basis
          ndim(i) = ncabs
          AO(i)   = AOdfCABS
       ELSE
          CALL LSQUIT('Error in get_AO_K in specification of AOs',-1)
       ENDIF
    ENDDO

    IF ((AO(3).NE.AORegular).OR.(AO(4).NE.AORegular)) &
         &  CALL LSQUIT('Error in get_AO_K in AO3 or AO4. Both should be AOregular',-1)

    IF (intspecConvert(5).EQ.'C') THEN
       Oper = coulombOperator
    ELSE
       CALL LSQUIT('Error in get_AO_K in specification of operator',-1)
    ENDIF

    ! Quick fix because we need to pass an array
    call mat_init(Karr(1),K%nrow,K%ncol)
    call mat_zero(Karr(1))
    doMPI = MyLsitem%SETTING%scheme%doMPI
    MyLsitem%SETTING%scheme%doMPI = .NOT.(((AO(1).EQ.AOdfCABS).OR.(AO(2).EQ.AOdfCABS)).OR.&
         & ((AO(3).EQ.AOdfCABS).OR.(AO(4).EQ.AOdfCABS)))
    CALL ii_get_exchange_mat_mixed(DECinfo%output,DECinfo%output,MyLsitem%SETTING,(/D/),ndmat,Dsym,&
         &                             Karr,AO(1),AO(3),AO(2),AO(4),Oper)
    MyLsitem%SETTING%scheme%doMPI = doMPI
    call mat_assign(K,Karr(1)) 
    call mat_free(Karr(1))

  end subroutine get_AO_K


#ifdef MOD_UNRELEASED
  !> Purpose: calculate AO int. in batches and transform them to
  !           full MO basis (non T1-transformed)
  !           The batches are then packed using permutational
  !           symmetry and are kept in memory (PDM if MPI)
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine get_t1_free_gmo(mo_ccsd,mylsitem,Co,Cv,pgmo_diag,pgmo_up, &
       & nb,no,nv,MOinfo)

    implicit none

    !> number of orbitals:
    integer, intent(in) :: nb, no, nv
    !> SCF transformation matrices:
    real(realk), pointer, intent(in) :: Co(:,:), Cv(:,:)
    !> performed MO-based CCSD calculation ?
    logical, intent(inout) :: mo_ccsd
    !> array with packed gmo on output:
    ! (intent in needed for the slaves)
    type(tensor), intent(inout) :: pgmo_diag, pgmo_up

    !> variables used for MO batch and integral transformation
    integer :: ntot ! total number of MO
    real(realk), pointer :: Cov(:,:), CP(:,:), CQ(:,:)
    real(realk), pointer :: w1(:), w2(:), gao(:)
    integer(kind=long)   :: w1_size, w2_size, gao_size
    integer :: Nbatch, PQ_batch, dimP, dimQ, idb, iub
    integer :: P_sta, P_end, Q_sta, Q_end
    type(MObatchInfo), intent(out) :: MOinfo
    logical :: local_moccsd

    !> variables used for AO batch construction and AO integral calculation
    integer :: alphaB, gammaB, dimAlpha, dimGamma
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
    integer :: iorb, idx, K
!ICHOR
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
!THERMITE
    integer, pointer :: batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchAlpha(:),orb2batchGamma(:)
    integer, pointer :: batchsizeGamma(:), batchindexGamma(:)
    ! Screening integrals stuff:
    type(DECscreenITEM) :: DecScreen
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    integer, pointer :: batchdimAlpha(:),batchdimGamma(:)
    Character :: INTSPEC(5)
    logical :: fullRHS, doscreen
    integer :: MaxAllowedDimAlpha, MaxActualDimAlpha, nbatchesAlpha
    integer :: MaxAllowedDimGamma, MaxActualDimGamma, nbatchesGamma

    !> CHECKING and MEASURING variables
    real(realk) :: tcpu, twall, time_start, timewall_start 
    logical :: print_debug

    ! MPI variables:
    logical :: master, local, gdi_lk, gup_lk, use_bg_buf
    integer(kind=8) :: nbu
    integer :: myload, win
    integer(kind=ls_mpik) :: ierr, myrank, nnod, dest
    integer, pointer      :: tasks(:)

    !> LS item with information needed for integrals
    type(lsitem), intent(inout) :: MyLsItem


    call time_start_phase(PHASE_WORK)
    ntot = no + nv
    use_bg_buf = mem_is_background_buf_init()
    if(use_bg_buf)then
       nbu = mem_get_bg_buf_free()
    else
       nbu = 0
    endif

    ! Set integral info
    ! *****************
    INTSPEC(1)  = 'R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)  = 'R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)  = 'R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)  = 'R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)  = 'C' !C = Coulomb operator
    IF(DECinfo%useIchor)THEN
       iprint = 0           !print level for Ichor Integral code
       MoTrans = .FALSE.    !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       iAO = 1
       call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
    ENDIF
    doscreen = MyLsItem%setting%scheme%cs_screen.OR. &
         & MyLsItem%setting%scheme%ps_screen
    ! Set MPI related info
    ! ********************
    master      = .true.
    local       = .true.
    myrank      = int(0,kind=ls_mpik)
    nnod        = 1
    !> logical stating if the windows of int. array are locked:
    gdi_lk      = .false.
    gup_lk      = .false.
#ifdef VAR_MPI
    myrank      = infpar%lg_mynum
    nnod        = infpar%lg_nodtot
    master      = (myrank == infpar%master)
#endif
    print_debug = (DECinfo%PL>2.or.DECinfo%cc_driver_debug.and.master)
    if (nnod>1) local=.false.

    ! Some timings
    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',time_start,timewall_start,DECinfo%output)

    ! Initialize stuff
    IF(DECinfo%useIchor)THEN
       nullify(AOGammabatchinfo)
       nullify(AOalphabatchinfo)    
    ELSE
       nullify(orb2batchAlpha)
       nullify(batchdimAlpha)
       nullify(batchsizeAlpha)
       nullify(batch2orbAlpha)
       nullify(batchindexAlpha)
       nullify(orb2batchGamma)
       nullify(batchdimGamma)
       nullify(batchsizeGamma)
       nullify(batch2orbGamma)
       nullify(batchindexGamma)
    ENDIF
    nullify(Cov)
    nullify(CP)
    nullify(CQ)
    nullify(w1)
    nullify(w2)
    nullify(MOinfo%DimInd1)
    nullify(MOinfo%DimInd2)
    nullify(MOinfo%StartInd1)
    nullify(MOinfo%StartInd2)
    nullify(tasks)

    !======================================================================
    !                      Get Dimension of batches                       !
    !======================================================================
    ! Get minimum mem. required in the MO-CCSD residual calculation
    if (master) then 

       call get_MO_and_AO_batches_size(mo_ccsd,local_moccsd,ntot,nb,no,nv, &
            & dimP,Nbatch,MaxAllowedDimAlpha,MaxAllowedDimGamma,MyLsItem,.false.)
       if (.not.mo_ccsd) return

       if (local_moccsd) then 
          write(DECinfo%output,*) 'Using MO-CCSD local scheme'
       else if (.not.local) then
          write(DECinfo%output,*) 'Using MO-CCSD: PDM scheme'
       else
          write(DECinfo%output,*) 'Using MO-CCSD: non-MPI scheme'
       end if
       if (print_debug) then
          write(DECinfo%output,'(a,I4,a,I4)') ' BATCH: Number of MO batches      = ', &
               & Nbatch*(Nbatch+1)/2, ' with maximum size', dimP
       end if

       ! Initialize gmo arrays:
       call init_gmo_arrays(ntot,dimP,Nbatch,local,local_moccsd,pgmo_diag,pgmo_up)

    end if
    !======================================================================



    !==================================================
    !                  Batch construction             !
    !==================================================

    call time_start_phase(PHASE_COMM)
    ! MPI: Waking slaves up:
#ifdef VAR_MPI
    StartUpSlaves: if(master.and.nnod>1) then
       call ls_mpibcast(CCGETGMO,infpar%master,infpar%lg_comm)
       call mpi_communicate_get_gmo_data(mo_ccsd,MyLsItem,Co,Cv, &
            & pgmo_diag,pgmo_up,nb,no,nv,Nbatch)
    endif StartUpSlaves

    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(dimP,infpar%master)
    call ls_mpi_buffer(Nbatch,infpar%master)
    call ls_mpi_buffer(local_moccsd,infpar%master)
    call ls_mpi_buffer(MaxAllowedDimAlpha,infpar%master)
    call ls_mpi_buffer(MaxAllowedDimGamma,infpar%master)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
#endif
    call time_start_phase(PHASE_WORK)


    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    IF(DECinfo%useIchor)THEN
       iAO = 4 !Gamma is the 4. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the MaxAllowedDimGamma, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
            & nbatchesGamma,DECinfo%output)
       call mem_alloc(AOGammabatchinfo,nbatchesGamma)
       !Construct the batches of AOS based on the MaxAllowedDimGamma, the requested
       !size of the AO batches - MaxAllowedDimGamma must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimGamma must be less og equal to MaxAllowedDimGamma
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
            & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
    ELSE
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchGamma,nb)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma, &
            & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma, &
            & nbatchesGamma,orb2BatchGamma,'R')
    ENDIF

    if (print_debug) write(DECinfo%output,'(a,I4,a,I4)') & 
         & ' BATCH: Number of Gamma batches   = ', nbatchesGamma, &
         & ' with maximum size', MaxActualDimGamma 


    IF(.NOT.DECinfo%useIchor)THEN
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbGamma,nbatchesGamma)
       do idx=1,nbatchesGamma
          call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx))
          batch2orbGamma(idx)%orbindex = 0
          batch2orbGamma(idx)%norbindex = 0
       end do
       do iorb=1,nb
          idx = orb2batchGamma(iorb)
          batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
          K = batch2orbGamma(idx)%norbindex
          batch2orbGamma(idx)%orbindex(K) = iorb
       end do
    endif


    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

    IF(DECinfo%useIchor)THEN
       iAO = 3 !Alpha is the 3. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the MaxAllowedDimAlpha, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
            & nbatchesAlpha,DECinfo%output)
       call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
       !Construct the batches of AOS based on the MaxAllowedDimAlpha, the requested
       !size of the AO batches - MaxAllowedDimAlpha must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimAlpha must be less og equal to MaxAllowedDimAlpha
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
            & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
    ELSE
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchAlpha,nb)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha, &
            & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha, &
            & nbatchesAlpha,orb2BatchAlpha,'R')
    ENDIF

    if (print_debug) write(DECinfo%output,'(a,I4,a,I4)') & 
         & ' BATCH: Number of Alpha batches   = ', nbatchesAlpha, &
         & ' with maximum size',MaxActualDimAlpha


    IF(.NOT.DECinfo%useIchor)THEN
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbAlpha,nbatchesAlpha)
       do idx=1,nbatchesAlpha
          call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
          batch2orbAlpha(idx)%orbindex = 0
          batch2orbAlpha(idx)%norbindex = 0
       end do
       do iorb=1,nb
          idx = orb2batchAlpha(iorb)
          batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
          K = batch2orbAlpha(idx)%norbindex
          batch2orbAlpha(idx)%orbindex(K) = iorb
       end do
    ENDIF

    ! **************************************
    ! * Allocate Memory to working arrays  *
    ! **************************************

    gao_size = int(i8*nb*nb*MaxActualDimAlpha*MaxActualDimGamma, kind=long)
    w1_size = max(ntot*MaxActualDimAlpha*MaxActualDimGamma, ntot*dimP*dimP)
    w1_size = int(i8*w1_size*ntot, kind=long)
    w2_size = max(ntot*nb*MaxActualDimAlpha*MaxActualDimGamma, &
       & MaxActualDimGamma*ntot*ntot*dimP, ntot*ntot*dimP*dimP)
    w2_size = int(i8*w2_size, kind=long)

    if(use_bg_buf)then
       if(w1_size+w2_size > nbu) then
          print *, "Warning(get_t1_free_gmo):  This should not happen, if the memory counting is correct"
       endif

       ! Get full MO coeficients:
       call mem_pseudo_alloc(Cov,i8*nb,i8*ntot)
       Cov(:,:no)       = Co
       Cov(:,no+1:ntot) = Cv
        
       ! CMO arrays:
       call get_MO_batches_info(MOinfo, dimP, ntot, Nbatch)
       call mem_pseudo_alloc(CP,i8*MaxActualDimAlpha,i8*dimP)
       call mem_pseudo_alloc(CQ,i8*MaxActualDimGamma,i8*dimP)

       call mem_pseudo_alloc(gao, gao_size)
       call mem_pseudo_alloc(w1, w1_size)
       call mem_pseudo_alloc(w2, w2_size)
    else
       ! Get full MO coeficients:
       call mem_alloc(Cov,nb,ntot)
       Cov(:,:no)       = Co
       Cov(:,no+1:ntot) = Cv
        
       ! CMO arrays:
       call get_MO_batches_info(MOinfo, dimP, ntot, Nbatch)
       call mem_alloc(CP,MaxActualDimAlpha,dimP)
       call mem_alloc(CQ,MaxActualDimGamma,dimP)

       call mem_alloc(gao, gao_size)
       call mem_alloc(w1, w1_size)
       call mem_alloc(w2, w2_size)
    endif


    ! Sanity checks for matrix sizes which need to be filled:
    if (max(w1_size, w2_size) > MaxInt) then
       call lsquit("ERROR(CCSD):matrix sizes too large, &
            & please recompile with 64bit integers",-1)
    endif

    IF(DECinfo%useIchor)THEN
       !Calculate Screening integrals 
       SameMOL = .TRUE. !Specifies same molecule on all centers 
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
    ELSE
       ! *******************************************************
       ! *  This subroutine builds the full screening matrix.
       call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting, &
            & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
       
       if (mylsitem%setting%scheme%cs_screen .or. mylsitem%setting%scheme%ps_screen) then
          call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
               & nb,nbatchesAlpha,nbatchesGamma,&
               & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
               & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
       end if
       ! *******************************************************
    ENDIF

#ifdef VAR_MPI
    ! Calculate the batches for a good load balance
    call mem_alloc(tasks,nbatchesAlpha*nbatchesGamma)

    myload = 0
    tasks  = 0
    IF(DECinfo%useIchor)THEN
       call mem_alloc(batchdimAlpha,nbatchesAlpha)
       do idx=1,nbatchesAlpha
          batchdimAlpha(idx) = AOAlphabatchinfo(idx)%dim 
       enddo
       call mem_alloc(batch2orbAlpha,nbatchesAlpha)
       do idx=1,nbatchesAlpha
          call mem_alloc(batch2orbAlpha(idx)%orbindex,1)
          batch2orbAlpha(idx)%orbindex(1) = AOAlphabatchinfo(idx)%orbstart
          batch2orbAlpha(idx)%norbindex = 1
       end do
       call mem_alloc(batchdimGamma,nbatchesGamma)
       do idx=1,nbatchesGamma
          batchdimGamma(idx) = AOGammabatchinfo(idx)%dim 
       enddo
       call mem_alloc(batch2orbGamma,nbatchesGamma)
       do idx=1,nbatchesGamma
          call mem_alloc(batch2orbGamma(idx)%orbindex,1)
          batch2orbGamma(idx)%orbindex(1) = AOGammabatchinfo(idx)%orbstart
          batch2orbGamma(idx)%norbindex = 1
       end do
       call distribute_mpi_jobs(tasks,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
            &batchdimGamma,myload,nnod,myrank,4,no,nv,nb,batch2orbAlpha,&
            &batch2orbGamma)
       call mem_dealloc(batchdimAlpha)
       call mem_dealloc(batchdimGamma)
       do idx=1,nbatchesAlpha
          call mem_dealloc(batch2orbAlpha(idx)%orbindex)
       end do
       call mem_dealloc(batch2orbAlpha)
       do idx=1,nbatchesGamma
          call mem_dealloc(batch2orbGamma(idx)%orbindex)
       end do
       call mem_dealloc(batch2orbGamma)
    ELSE
       call distribute_mpi_jobs(tasks,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
            &batchdimGamma,myload,nnod,myrank,4,no,nv,nb,batch2orbAlpha,&
            &batch2orbGamma)
    ENDIF
#endif
    myload = 0

    fullRHS = (nbatchesGamma.eq.1).and.(nbatchesAlpha.eq.1)

    !**********************************
    ! Begin the loop over gamma batches
    !**********************************


    BatchGamma: do gammaB = 1,nbatchesGamma            ! batches of AO batches
       IF(DECinfo%useIchor)THEN
          dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
          GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
          GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
          AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch needed by MAIN_ICHORERI_DRIVER
          AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch needed by MAIN_ICHORERI_DRIVER
       ELSE
          dimGamma   = batchdimGamma(gammaB)                         ! Dimension of gamma batch
          GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
          GammaEnd   = batch2orbGamma(gammaB)%orbindex(dimGamma)     ! Last index in gamma batch
       ENDIF
       
       !**********************************
       ! Begin the loop over alpha batches
       !**********************************

       BatchAlpha: do alphaB = 1,nbatchesAlpha         ! batches of AO batches

          if (nnod>1) then 
             ! check if the current job is to be done by current node
             if (tasks(alphaB + (gammaB-1)*nbatchesAlpha)/=myrank) cycle
          end if

          IF(DECinfo%useIchor)THEN
             dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
             AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
             AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
             AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
             AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
          ELSE
             dimAlpha   = batchdimAlpha(alphaB)                        ! Dimension of alpha batch
             AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)           ! First index in alpha batch
             AlphaEnd   = batch2orbAlpha(alphaB)%orbindex(dimAlpha)    ! Last index in alpha batch
          ENDIF
          myload     = myload + dimAlpha*dimGamma


          IF(DECinfo%useIchor)THEN
             call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nb,nb,dimAlpha,dimGamma,&
                  & gao,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
                  & AOGammaStart,AOGammaEnd,MoTrans,nb,nb,dimAlpha,dimGamma,NoSymmetry,DECinfo%IntegralThreshold)
          ELSE
             ! setup RHS screening - here we only have a set of AO basisfunctions
             !                      so we use the batchscreening matrices.
             !                      like BatchfilenamesCS(alphaB,gammaB)
             ! Note that it is faster to calculate the integrals in the form
             ! (dimAlpha,dimGamma,nbasis,nbasis) so the subset of the AO basis is used on the LHS
             ! but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)
             IF(doscreen) Mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
             IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p
             
             ! Get AO integrals using (beta,delta,alphaB,gammaB) ordering
             ! **********************************************************
             call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
                  & Mylsitem%setting,gao,batchindexAlpha(alphaB),batchindexGamma(gammaB), &
                  & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nb,nb,dimAlpha, &
                  & dimGamma,fullRHS,INTSPEC,DECinfo%IntegralThreshold)
          ENDIF
          idb = 0
          iub = 0
          ! Loop over MO batches:
          BatchPQ: do PQ_batch = 1, MOinfo%nbatch


             P_sta  = MOinfo%StartInd1(PQ_batch)
             dimP   = MOinfo%DimInd1(PQ_batch)
             Q_sta  = MOinfo%StartInd2(PQ_batch)
             dimQ   = MOinfo%DimInd2(PQ_batch)

             ! w2 contains MO integral batch on output
             call gao_to_gmo(gao,w1,w2,Cov,CP,CQ,nb,ntot,AlphaStart,dimAlpha, &
                  & GammaStart,dimGamma,P_sta,dimP,Q_sta,dimQ,pgmo_diag, &
                  & pgmo_up,gdi_lk,gup_lk,win,dest)

             if (P_sta==Q_sta) then
                idb = idb + 1 
                if (.not.local) then
                   !LOCK WINDOW AND LOCK_SET = .true.
                   win = idb
#ifdef VAR_MPI
                   call tensor_lock_win(pgmo_diag,win,'s')
                   gdi_lk = .true. 
#endif
                end if
                call pack_and_add_gmo(w2,pgmo_diag,idb,ntot,dimP,dimQ,.true.,w1)
             else 
                iub = iub + 1 
                if (.not.local) then
                   !LOCK WINDOW AND LOCK_SET = .true.
                   win = iub
#ifdef VAR_MPI
                   call tensor_lock_win(pgmo_up,win,'s')
                   gup_lk = .true.
#endif
                end if
                call pack_and_add_gmo(w2,pgmo_up,iub,ntot,dimP,dimQ,.false.,w1)
             end if

          end do BatchPQ


       end do BatchAlpha
    end do BatchGamma

    ! Free integral stuff
    ! *******************
    IF(DECinfo%useIchor)THEN
       call FREE_SCREEN_ICHORERI()
       call mem_dealloc(AOGammabatchinfo)
       call mem_dealloc(AOAlphabatchinfo)
    ELSE
       nullify(Mylsitem%setting%LST_GAB_LHS)
       nullify(Mylsitem%setting%LST_GAB_RHS)
       call free_decscreen(DECSCREEN)
       
       ! Free gamma stuff
       call mem_dealloc(orb2batchGamma)
       call mem_dealloc(batchdimGamma)
       call mem_dealloc(batchsizeGamma)
       call mem_dealloc(batchindexGamma)
       do idx=1,nbatchesGamma
          call mem_dealloc(batch2orbGamma(idx)%orbindex)
          batch2orbGamma(idx)%orbindex => null()
       end do
       call mem_dealloc(batch2orbGamma)
       
       ! Free alpha stuff
       call mem_dealloc(orb2batchAlpha)
       call mem_dealloc(batchdimAlpha)
       call mem_dealloc(batchsizeAlpha)
       call mem_dealloc(batchindexAlpha)
       do idx=1,nbatchesAlpha
          call mem_dealloc(batch2orbAlpha(idx)%orbindex)
          batch2orbAlpha(idx)%orbindex => null()
       end do
       call mem_dealloc(batch2orbAlpha)
    ENDIF

#ifdef VAR_MPI
    ! UNLOCK REMAINING WINDOWS
    if (gdi_lk) then
       call tensor_unlock_win(pgmo_diag,win)
    else if (gup_lk) then
       call tensor_unlock_win(pgmo_up,win)
    end if

    call mem_dealloc(tasks)
    ! Problem specific to one sided comm. and maybe bcast,
    ! We must use a barrier after one sided communication epoc:
    if (.not.local_moccsd) then
       call time_start_phase(PHASE_IDLE)
       call lsmpi_barrier(infpar%lg_comm)
       call time_start_phase(PHASE_WORK)
    end if
#endif

    ! Free matrices:
    if(use_bg_buf)then
       call mem_pseudo_dealloc(w2)
       call mem_pseudo_dealloc(w1)
       call mem_pseudo_dealloc(gao)
       call mem_pseudo_dealloc(CQ)
       call mem_pseudo_dealloc(CP)
       call mem_pseudo_dealloc(Cov)
    else
       call mem_dealloc(w2)
       call mem_dealloc(w1)
       call mem_dealloc(gao)
       call mem_dealloc(CQ)
       call mem_dealloc(CP)
       call mem_dealloc(Cov)
    endif

    call LSTIMER('get_t1_free_gmo',tcpu,twall,DECinfo%output)

  end subroutine get_t1_free_gmo


  !> Purpose: Calculate MO batches size based on available memory
  !           and memory requirements in MO-CCSD residual routine.
  !           Calculate Max. AO batches based on MO batches size,
  !           Min. AO batches size, available memory and mem. 
  !           requirements in get_t1_free_gmo routine.
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  subroutine get_MO_and_AO_batches_size(mo_ccsd,local,ntot,nb,no,nv, &
       & dimMO,Nbatch,MaxAlpha,MaxGamma,MyLsItem,mpi_split)

    implicit none

    !> performed MO-based CCSD calculation ?
    logical, intent(inout) :: mo_ccsd, local
    !> number of orbitals:
    integer, intent(in) :: ntot, nb, no, nv
    !> MO batches stuff:
    integer, intent (inout) :: dimMO, Nbatch
    !> AO batches stuff:
    integer, intent (inout) :: MaxAlpha, MaxGamma
    type(lsitem), intent(inout) :: MyLsItem
    logical, intent(in) :: mpi_split

    real(realk) :: MemNeed, MemFree, factor
    integer(kind=long) :: min_mem
    integer :: MinAOBatch, MinMOBatch, na, ng, nnod, magic,iAO

    MinMOBatch = min(15,ntot)
    factor = 0.8E0_realk
    dimMO = MinMOBatch
    local = .false.
    nnod  = 1
#ifdef VAR_MPI
    nnod  = infpar%lg_nodtot
#endif

    !===========================================================
    ! Get MO batch size depending on MO-ccsd residual routine.
    call get_currently_available_memory(MemFree)

    if(mem_is_background_buf_init())then
       MemFree = factor*MemFree + (dble(mem_get_bg_buf_n())*8.0E0_realk)/(1024.0**3)
    endif

    if (nnod>1) then

       ! SELECT SCHEME (storage of MO int.): 
       !
       ! 0-4 are reserved to standard CCSD (Patrick's code)
       ! 
       ! 5: Local scheme: more memory required but no one 
       !    sided communication.
       !
       ! 6: PDM scheme: batches are distributed in PDM using
       !    one sided communication.

       if (DECinfo%force_scheme.and.DECinfo%en_mem==6) then
          print *,"!!FORCING MO-CCSD LOCAL SCHEME!!"
          local = .true.
       else if (DECinfo%force_scheme.and.DECinfo%en_mem==5) then
          print *,"!!FORCING MO-CCSD PDM SCHEME!!"
          local = .false.
       else
          ! try first for scheme with highest requirements --> fastest
          local = .true.
          MemNeed = get_mem_MO_CCSD_residual(local,nb,no,nv,dimMO)

          ! if not enough mem then switch to full PDM scheme:
          if (MeMNeed>factor*MemFree) then
             local = .false.
          end if
       end if
    end if

    MemNeed = get_mem_MO_CCSD_residual(local,nb,no,nv,dimMO)

    do while ((MemNeed<factor*MemFree).and.(dimMO<=ntot))
       dimMO = dimMO + 1
       MemNeed = get_mem_MO_CCSD_residual(local,nb,no,nv,dimMO)
    end do

    if (dimMO>=ntot) then
       dimMO = ntot
    else if (dimMO<=MinMOBatch) then
       dimMO = MinMOBatch
    else
       dimMO = dimMO - 1
    end if

    ! mpi_split should be true when we want to estimate the workload associated
    ! to a DEC fragment and eventually split the slots. In this case, the next
    ! step must be skiped.
    if (.not.mpi_split) then
       ! Check that every nodes will have a job in residual calc.
       ! But the dimension of the batch must stay above MinMOBatch.
       magic  = int(1.5*nnod)
       Nbatch = ((ntot-1)/dimMO+1)
       Nbatch = Nbatch*(Nbatch+1)/2

       do while (Nbatch<magic.and.(dimMO>MinMOBatch).and.nnod>1)
          dimMO = dimMO-1
          Nbatch = ((ntot-1)/dimMO+1)
          Nbatch = Nbatch*(Nbatch+1)/2
          if (dimMO<MinMOBatch) then
             dimMO = MinMOBatch
             exit
          end if
       end do
    end if
    ! sanity check:
    MemNeed = get_mem_MO_CCSD_residual(local,nb,no,nv,dimMO)
    if ((MemFree-MemNeed)<=0.0E0_realk) then
       mo_ccsd = .false.
       write(DECinfo%output,'(a,F12.5,a)') '   Available memory:',MemFree,' GB'
       write(DECinfo%output,'(a,F12.5,a)') '   Required memory :',MemNeed,' GB'
       if (DECinfo%force_scheme) then
          call lsquit('Insufficient memory in MO-based CCSD (remove force scheme)',DECinfo%output)
       else
          write(DECinfo%output,*) 'WARNING: Insufficient memory in MO-based CCSD, &
               & back to standard algorithm.'
       end if
       return
    end if
    Nbatch = (ntot-1)/dimMO + 1

    !===========================================================
    ! Get AO batche size depending on get_t1_free_gmo routine.
    IF(DECinfo%useIchor)THEN
       !Determine the minimum allowed AObatch size MinAObatch
       !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
       !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
       !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
       !'R'  !Specifies that it is the Regular AO basis that should be batched
       iAO = 1 !the center that the batching should occur on.  
       call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
    ELSE
       call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
    ENDIF

    MemNeed = get_mem_t1_free_gmo(local,nb,no,nv,dimMO,MinAObatch,MinAObatch)

    MaxGamma = MinAObatch
    MaxAlpha = MinAObatch
    do while ((MemNeed<factor*MemFree).and.(MaxGamma<=nb)) 
       MaxGamma = MaxGamma + 1
       MemNeed = get_mem_t1_free_gmo(local,nb,no,nv,dimMO,MaxAlpha,MaxGamma)
    end do
    if (MaxGamma>=nb) then
       MaxGamma = nb
    else if (MaxGamma<=MinAObatch) then
       MaxGamma = MinAObatch
    else 
       MaxGamma = MaxGamma - 1
    end if
    do while ((MemNeed<factor*MemFree).and.(MaxAlpha<=nb)) 
       MaxAlpha = MaxAlpha + 1
       MemNeed = get_mem_t1_free_gmo(local,nb,no,nv,dimMO,MaxAlpha,MaxGamma)
    end do
    if (MaxAlpha>=nb) then
       MaxAlpha = nb
    else if (MaxAlpha<=MinAObatch) then
       MaxAlpha = MinAObatch
    else 
       MaxAlpha = MaxAlpha - 1
    end if

    ! mpi_split should be true when we want to estimate the workload associated
    ! to a DEC fragment and eventually split the slots. In this case, the next
    ! step must be skiped.
    if (.not.mpi_split) then
       ! Check that every nodes has a job:
       magic = int(2*nnod)
       ng    = ((nb-1)/MaxGamma+1)
       na    = ((nb-1)/MaxAlpha+1)

       ! Number of Alpha batches must be at least magic
       if (na*ng<magic.and.(MaxAlpha>MinAObatch).and.nnod>1)then
          MaxAlpha = (nb/magic)
          if (MaxAlpha<MinAObatch) MaxAlpha = MinAObatch
       end if

       na    = ((nb-1)/MaxAlpha+1)
       if (na*ng<magic.and.(MaxAlpha==MinAObatch).and.nnod>1)then
          do while(na*ng<magic)
             MaxGamma = MaxGamma - 1
             if (MaxGamma<MinAObatch) then
                MaxGamma = MinAObatch
                exit
             end if
             ng    = ((nb-1)/MaxGamma+1)
          end do
       endif
    end if

    ! sanity check:
    MemNeed = get_mem_t1_free_gmo(local,nb,no,nv,dimMO,MaxAlpha,MaxGamma)
    if ((MemFree-MemNeed)<=0.0E0_realk) then
       mo_ccsd = .false.
       write(DECinfo%output,*) 'WARNING: Insufficient memory in MO-based CCSD, &
            & back to standard algorithm.'
       write(DECinfo%output,'(a,F12.5,a)') '   Available memory:',MemFree,' GB'
       write(DECinfo%output,'(a,F12.5,a)') '   Required memory :',MemNeed,' GB'
       return
    end if

    ! Check that there is enough mem in bg buffer:
    if(mem_is_background_buf_init())then
       MemNeed = get_mem_t1_free_gmo(local,nb,no,nv,dimMO,MaxAlpha,MaxGamma,2)
       MemNeed = MemNeed + get_mem_MO_CCSD_residual(local,nb,no,nv,dimMO,2)
       MemFree = (dble(mem_get_bg_buf_n())*8.0E0_realk)/(1024.0**3)
       if ((MemFree-MemNeed)<=0.0E0_realk) then
          mo_ccsd = .false.
          write(DECinfo%output,*) 'WARNING: Insufficient memory in background buffer for &
             & MO-based CCSD, back to standard algorithm.'
          write(DECinfo%output,'(a,F12.5,a)') '   Available memory:',MemFree,' GB'
          write(DECinfo%output,'(a,F12.5,a)') '   Required memory :',MemNeed,' GB'
          return
       end if
    endif

  end subroutine get_MO_and_AO_batches_size


  !> Purpose: Get memory required in get_t1_free_gmo depending on AO 
  !           batch dimension.
  !           Get min. required memory: AlphaDim = GammaDim = MinDimAO
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  function get_mem_t1_free_gmo(local,N,O,V,X,adim,gdim,which_mem) result(mem_out)

    implicit none 

    ! N: tot number of AO
    ! O: number of occ. orbs.
    ! V: number of virt. orbs.
    ! X: dimension of MO batch.
    integer,  intent(in) :: N, O, V, X
    !> use local scheme?
    logical :: local
    !> AO stuff:
    integer, intent(in) :: adim, gdim
    !> memory needed:
    real(realk) :: mem_out
    !> define memory output:
    !  which_mem = 1 (default) returns total memory
    !  which_mem = 2 returns bg memory only
    !  which_mem = 3 returns non bg memory only
    integer, intent(in), optional :: which_mem
    !> intermediate memory:
    integer :: nnod, nMOB, M, wm
    integer(kind=long) :: mem0, mem1, mem2
    integer(kind=long) :: ntile_max, mem_bg, mem_not_bg

    wm = 1
    if (present(which_mem)) wm = which_mem

    nnod = 1
    M = O+V
    nMOB = (M-1)/X + 1
    mem_out = 0.0_realk
    mem_bg = 0.0_realk
    mem_not_bg = 0.0_realk

#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
#endif
    if (local) nnod = 1

    ! Transfo. matrices:
    mem_bg = mem_bg + N*M + adim*X + gdim*X

    ! AO stuff:
    mem_not_bg = mem_not_bg + 4*N + N*N*adim*gdim 

    ! Packed gmo diag blocks:
    ntile_max = (nMOB-1)/nnod + 3
    mem_bg = mem_bg + ntile_max*X*(X+1)*M*(M+1)/4
    ! Packed gmo upper blocks:
    ntile_max = (nMOB*(nMOB-1)/2 - 1)/nnod + 3
    mem_bg = mem_bg + ntile_max*X*X*M*(M+1)/2

    ! MO stuff:
    mem_not_bg = mem_not_bg + 5*nMOB*nMOB + 1

    ! Working arrays:
    mem_bg = mem_bg + max(M*N*adim*gdim, M*M*gdim*X)
    mem_bg = mem_bg + max(M*M*adim*gdim, M*M*X*X)

    select case(wm)
    case(1)
      mem_out = (mem_bg + mem_not_bg)*8.0E0_realk/(1.024E3_realk**3) 
    case(2)
      mem_out = (mem_bg)*8.0E0_realk/(1.024E3_realk**3) 
    case(3)
      mem_out = (mem_not_bg)*8.0E0_realk/(1.024E3_realk**3) 
    case default
       call lsquit("ERROR(get_mem_t1_free_gmo): input not valid",DECinfo%output)
    end select

  end function get_mem_t1_free_gmo


  !> Purpose: Get memory required in get_ccsd_residual_mo_ccsd 
  !           depending on MO batch dimension.
  !           Get min. required memory when X is set to 1 
  !           In total this should amount to ~ 9*v2o2
  !           + M**4 (scheme 6) or M**4/nnodes (scheme 5)
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  function get_mem_MO_CCSD_residual(local,N,O,V,X,which_mem) result(mem_out)

    implicit none 

    ! N: tot number of AO
    ! O: number of occ. orbs.
    ! V: number of virt. orbs.
    ! X: dimension of MO batch.
    integer,  intent(in) :: N, O, V, X
    !> use local scheme?
    logical, intent(in) :: local
    !> memory needed:
    real(realk) :: mem_out
    !> define memory output:
    !  which_mem = 1 (default) returns total memory
    !  which_mem = 2 returns bg memory only
    !  which_mem = 3 returns non bg memory only
    integer, intent(in), optional :: which_mem
    !> intermediate memory:
    integer :: nnod, nMOB, M, wm
    integer(kind=long) :: mem0, mem1, mem2
    integer(kind=long) :: ntile_max, mem_bg, mem_not_bg

    wm = 1
    if (present(which_mem)) wm = which_mem

    M = O+V
    nMOB = (M-1)/X + 1
    nnod = 1
    mem_out = 0.0_realk
    mem_bg = 0.0_realk
    mem_not_bg = 0.0_realk
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
#endif
    if (local) nnod = 1

    ! Tiled array are kept as dense during the residual calc:
    ! (govov amplitudes and residual)
    mem_bg = mem_bg + 3*O*O*V*V

    ! Packed gmo diag blocks:
    ntile_max = (nMOB-1)/nnod + 3
    mem_bg = mem_bg + ntile_max*X*(X+1)*M*(M+1)/4

    ! Packed gmo upper blocks:
    ntile_max = (nMOB*(nMOB-1)/2 - 1)/nnod + 3
    mem_bg = mem_bg + ntile_max*X*X*M*(M+1)/2

    ! Working arrays (~ 3*v2o2):
    call get_mem_mo_ccsd_warrays(o,v,x, mem0, mem1, mem2)
    mem_bg = mem_bg + mem0 + mem1 + mem2

    ! T1-Transfo. matrices:
    mem_bg = mem_bg + V*M + O*M

    ! Batch of MO int:
    mem_bg = mem_bg + X*X*M*M

    ! Intermediates (B2prep, u2, G_Pi, H_aQ):
    mem_bg = mem_bg + O**4 + O*O*V*V + M*O + V*M

    ! T1-transformed integrals:
    mem_bg = mem_bg + 2*O*O*V*V 

    select case(wm)
    case(1)
      mem_out = (mem_bg + mem_not_bg)*8.0E0_realk/(1.024E3_realk**3) 
    case(2)
      mem_out = (mem_bg)*8.0E0_realk/(1.024E3_realk**3) 
    case(3)
      mem_out = (mem_not_bg)*8.0E0_realk/(1.024E3_realk**3) 
    case default
       call lsquit("ERROR(get_mem_MO_CCSD_residual): input not valid",DECinfo%output)
    end select

  end function get_mem_MO_CCSD_residual


  !> Purpose: Get maximum used memory for working array in the MO-based 
  !           CCSD residual calculation.
  !
  !> Author:  Pablo Baudin
  !> Date:    July 2015
  subroutine get_mem_mo_ccsd_warrays(o,v,x, mem0, mem1, mem2)

    implicit none 

    ! O: number of occ. orbs.
    ! V: number of virt. orbs.
    ! X: dimension of MO batch.
    integer,  intent(in) :: O, V, X
    !> Max memory used for working arrays:
    integer(kind=long), intent(out) :: mem0, mem1, mem2
    !> intermediates:
    integer:: m, vbar, obar

    m = v+o
    vbar = min(x,v)
    obar = min(x,o)

    ! I know many terms can be easily removed but it's good to keep 
    ! a one to one  mapping between the code and this routine (and the notes).

    ! Max memory used for tmp0:
    mem0 = max(x*m*x*m, v*vbar*o*o, v*v*obar*o, o*o*v*x, x*vbar*o*m, &
       & v*vbar*o*o, v*v*o*o)
    mem0 = int(i8*mem0, kind=long)

    ! Max memory used for tmp1:
    mem1 = max(x*x*m*o, m*o*v*o, m*obar*o*x, obar*o*o*o, vbar*v*o*o/2, &
       & obar*o*o*o/2, m*o*o*x/2, v*v*o*o/2, x*v*vbar*o, v*obar*o*x, &
       & x*o*v*x, x*vbar*o*o, v*v*o*o)
    mem1 = int(i8*mem1, kind=long)

    ! Max memory used for tmp2:
    mem2 = max(x*m*o*v, o*v*o*v, m*obar*o*o, x*m*o*o/2, v*o*o*x/2, v*v*o*o)
    mem2 = int(i8*mem2, kind=long)

  end subroutine get_mem_mo_ccsd_warrays


  !> Purpose: Initialization of arrays for MO integrals: 
  !           if NO MPI then the arrays are standard
  !           if MPI and local_moccsd then the arrays are RTAR
  !           i.e. all the tiles are stored on all the nodes
  !           if MPI and not local_moccsd then the arrays are TDAR
  !           i.e. the tiles are distributed among the nodes.
  !
  !> Author:  Pablo Baudin
  !> Date:    February 2014
  subroutine init_gmo_arrays(ntot,dimMO,Nbat,mpi,local_moccsd,pgmo_diag,pgmo_up)

    implicit none

    !> dimension parameters: 
    integer, intent(in) :: ntot, dimMO, Nbat
    !> logical for the type of arrays:
    logical, intent(in) :: mpi, local_moccsd
    !> gmo arrays:
    type(tensor), intent(inout) :: pgmo_diag, pgmo_up

    logical :: use_bg_buf
    character(4) :: at
    integer :: pgmo_dims

    use_bg_buf = mem_is_background_buf_init()

    ! define type of array:
    if (local_moccsd) then
       at = 'RTAR'
    else 
       at = 'TDAR'
    end if

    ! Declare one array for the diagonal batches:
    pgmo_dims = ntot*(ntot+1)*dimMO*(dimMO+1)/4
    call tensor_minit(pgmo_diag,[pgmo_dims,Nbat],2,local=mpi,atype=at,tdims=[pgmo_dims,1], &
       & bg=use_bg_buf)
    call tensor_zero(pgmo_diag)

    ! Declare one array for the upper diagonal batches if necesarry
    if (Nbat>1) then
      pgmo_dims = ntot*(ntot+1)*dimMO*dimMO/2
      call tensor_minit(pgmo_up,[pgmo_dims,Nbat*(Nbat-1)/2],2,local=mpi,atype=at, &
         & tdims=[pgmo_dims,1], bg=use_bg_buf)
      call tensor_zero(pgmo_up)
    end if

  end subroutine init_gmo_arrays


  !> Purpose: Transform AO int. into MO int. in batches
  !           
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine gao_to_gmo(gao,w1,w2,Cov,CP,CQ,nb,ntot,AlphaStart,dimAlpha, &
       & GammaStart,dimGamma,P_sta,dimP,Q_sta,dimQ,pgmo_diag,pgmo_up, &
       & gdi_lk,gup_lk,win,dest)

    implicit none

    integer, intent(in) :: nb, AlphaStart, dimAlpha, GammaStart, dimGamma
    integer, intent(in) :: ntot, P_sta, dimP, Q_sta, dimQ
    !> AO integral batch:
    real(realk), intent(inout) :: gao(:)
    !> working arrays, w2 contains MO integral batch on output:
    real(realk), intent(inout) :: w1(:), w2(:)
    real(realk), intent(in) :: Cov(nb,ntot)
    !> MPI related:
    type(tensor), intent(in) :: pgmo_diag, pgmo_up
    logical, intent(inout)  :: gdi_lk, gup_lk
    integer, intent(in) :: win
    integer(kind=ls_mpik), intent(in) :: dest

    real(realk) :: CP(dimAlpha,dimP), CQ(dimGamma,dimQ)
    integer :: AlphaEnd, GammaEnd, P_end, Q_end

    AlphaEnd = AlphaStart+dimAlpha-1
    GammaEnd = GammaStart+dimGamma-1
    P_end     = P_sta+dimP-1
    Q_end     = Q_sta+dimQ-1

    ! initialisation of transfo. matrices:
    CP = Cov(AlphaStart:AlphaEnd,P_sta:P_end)
    CQ = Cov(GammaStart:GammaEnd,Q_sta:Q_end)

    ! transfo Beta to r => [delta alphaB gammaB, r]
    call dgemm('t','n',nb*dimAlpha*dimGamma,ntot,nb,1.0E0_realk, &
         & gao,nb,Cov,nb,0.0E0_realk,w2,nb*dimAlpha*dimGamma)

#ifdef VAR_MPI
    ! UNLOCK WINDOW IF (LOCK_SET)
    if (gdi_lk) then
       call tensor_unlock_win(pgmo_diag,win)
       gdi_lk = .false.
    else if (gup_lk) then
       call tensor_unlock_win(pgmo_up,win)
       gup_lk = .false.
    end if
#endif

    ! transfo delta to s => [alphaB gammaB r, s]
    call dgemm('t','n',dimAlpha*dimGamma*ntot,ntot,nb,1.0E0_realk, &
         & w2,nb,Cov,nb,0.0E0_realk,w1,dimAlpha*dimGamma*ntot)

    ! transfo alphaB to P_batch => [gammaB r s, P]
    call dgemm('t','n',dimGamma*ntot*ntot,dimP,dimAlpha,1.0E0_realk, &
         & w1,dimAlpha,CP,dimAlpha,0.0E0_realk,w2,dimGamma*ntot*ntot)

    ! transfo gammaB to Q_batch => [r s P, Q]
    call dgemm('t','n',ntot*ntot*dimP,dimQ,dimGamma,1.0E0_realk, &
         & w2,dimGamma,CQ,dimGamma,0.0E0_realk,w1,ntot*ntot*dimP)

    ! transpose matrix => [P_batch, Q_batch, r, s]
    call mat_transpose(ntot*ntot,dimP*dimQ,1.0E0_realk,w1,0.0E0_realk,w2)

  end subroutine gao_to_gmo


  !> Purpose: Get information regarding MO batches for two 
  !           equivalent nested loops:
  !
  !> Author:  Pablo Baudin
  !> Date:    November 2013
  subroutine get_MO_batches_info(PQbatchInfo, dimBat, TotSize, Nbat)

    implicit none

    integer, intent(in) :: TotSize, dimBat, Nbat
    type(MObatchInfo), intent(inout) :: PQbatchInfo

    integer :: PQ_batch, Pbatch, Qbatch, dimP,  P_sta, dimQ, Q_sta, Njob, idb, iub

    ! short cuts:
    PQbatchInfo%nbatch = Nbat*(Nbat+1)/2
    Njob = PQbatchInfo%nbatch

    ! Allocate arrays:
    call mem_alloc(PQbatchInfo%StartInd1, Njob)
    call mem_alloc(PQbatchInfo%StartInd2, Njob)
    call mem_alloc(PQbatchInfo%dimInd1,   Njob)
    call mem_alloc(PQbatchInfo%dimInd2,   Njob)
    call mem_alloc(PQbatchInfo%dimTot,    Njob)
    call mem_alloc(PQbatchInfo%tileInd, Njob,2)

    ! Initialization
    PQ_batch = 1
    Q_sta = 1
    dimP = dimBat
    dimQ = dimBat
    idb = 0
    iub = 0

    ! Loop over MO batches:
    BatchQ: do Qbatch = 1, Nbat

       ! get dimension of last Q batch:
       if (Qbatch == Nbat) dimQ = TotSize - Q_sta + 1

       P_sta = 1
       BatchP: do Pbatch = 1, Qbatch

          ! get dimension of last P batch:
          if (Pbatch == Nbat) dimP = TotSize - P_sta + 1

          ! Store info about this PQ batch:
          PQbatchInfo%StartInd1(PQ_batch) = P_sta
          PQbatchInfo%StartInd2(PQ_batch) = Q_sta
          PQbatchInfo%dimInd1(PQ_batch)   = dimP
          PQbatchInfo%dimInd2(PQ_batch)   = dimQ

          ! DimTot contains the total dimension 
          if (P_sta==Q_sta) then
             idb = idb + 1
             PQbatchInfo%tileInd(PQ_batch,1) = idb
             PQbatchInfo%tileInd(PQ_batch,2) = 0
             PQbatchInfo%dimTot(PQ_batch) = dimP*dimQ
          else 
             iub = iub + 1
             PQbatchInfo%tileInd(PQ_batch,1) = iub
             PQbatchInfo%tileInd(PQ_batch,2) = 1
             PQbatchInfo%dimTot(PQ_batch) = 2*dimP*dimQ
          end if

          PQ_batch = PQ_batch + 1
          P_sta = P_sta + dimP

       end do BatchP
       ! restore dimension of P batch
       dimP = dimQ
       Q_sta = Q_sta + dimQ
    end do BatchQ

  end subroutine get_MO_batches_info


  !> Purpose: Pack MO integrals using symmetry of charge distribution.
  !           The contributions are summing over AO batch
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  subroutine pack_and_add_gmo(w2,pack_gmo,tile,ntot,dimP,dimQ,diag,w1)

    implicit none

    !> array with one batch of partial MO int.:
    real(realk), intent(in) :: w2(:)
    !> array containing the previous contributions
    !  to this MO int. batch, packed.
    type(tensor), intent(inout) :: pack_gmo
    !> index corresponding to the current batch:
    integer, intent(in) :: tile
    !> dimensions of array:
    integer, intent(in) :: ntot, dimP, dimQ
    !> Diagonal block ?
    logical, intent(in) :: diag
    !> working array:
    real(realk), intent(inout) :: w1(:)

    integer :: s, r, rs, q, ibatch, ipack, nnod
    integer(kind=long) :: ncopy

    nnod = 1
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
#endif

    if (diag) then
       ! 1st case: current batch corresponds to diagonal block, we 
       !           keep only the upper triangular part of the batch.

       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ntot,w2,w1,dimP,dimQ)&
       !$OMP PRIVATE(q,r,s,rs,ibatch,ipack)
       do s=1,ntot
          do r=1,s
             rs = r + (s-1)*ntot
             do q=1,dimQ
                ibatch = 1 + (q-1)*dimP + (rs-1)*dimP*dimQ
                ipack = 1 + q*(q-1)/2 + (s*(s-1)/2 + r-1)*(dimQ*(dimQ+1)/2)
                call dcopy(q,w2(ibatch),1,w1(ipack),1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       ncopy = dimQ*(dimQ+1)/2 + (ntot*(ntot-1)/2 + ntot-1)*(dimQ*(dimQ+1)/2)
       ! accumulate tile
       if (nnod>1.and.pack_gmo%itype==TT_TILED_DIST) then
          call time_start_phase(PHASE_COMM)
          call tensor_accumulate_tile(pack_gmo,tile,w1(1:ncopy),ncopy,lock_set=.true.)
          call time_start_phase(PHASE_WORK)
       else if (nnod>1.and.pack_gmo%itype==TT_TILED_REPL) then
          call daxpy(ncopy,1.0E0_realk,w1,1,pack_gmo%ti(tile)%t(:),1)
       else
          call daxpy(ncopy,1.0E0_realk,w1,1,pack_gmo%elm2(:,tile),1)
       end if

    else
       ! 2nd case: current batch corresponds to an upper diagonal block,
       !           we keep all the pq part and reduced r<=s.

       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ntot,w2,w1,dimP,dimQ)&
       !$OMP PRIVATE(r,s,rs,ibatch,ipack,ncopy)
       do s=1,ntot
          do r=1,s
             rs = r + (s-1)*ntot
             ibatch = 1 + (rs-1)*dimP*dimQ
             ipack = 1 + (s*(s-1)/2 + r-1)*dimP*dimQ
             ncopy = dimP*dimQ
             call dcopy(ncopy,w2(ibatch),1,w1(ipack),1)
          end do
       end do
       !$OMP END PARALLEL DO

       ncopy = (ntot*(ntot-1)/2 + ntot)*dimP*dimQ
       ! accumulate tile
       if (nnod>1.and.pack_gmo%itype==TT_TILED_DIST) then
          call time_start_phase(PHASE_COMM)
          call tensor_accumulate_tile(pack_gmo,tile,w1(1:ncopy),ncopy,lock_set=.true.)
          call time_start_phase(PHASE_WORK)
       else if (nnod>1.and.pack_gmo%itype==TT_TILED_REPL) then
          call daxpy(ncopy,1.0E0_realk,w1,1,pack_gmo%ti(tile)%t(:),1)
       else
          call daxpy(ncopy,1.0E0_realk,w1,1,pack_gmo%elm2(:,tile),1)
       end if

    end if

  end subroutine pack_and_add_gmo


  !> Purpose: Unpack MO integrals using symmetry of charge distribution.
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  subroutine unpack_gmo(gmo,pack_gmo,tile,ntot,dimP,dimQ,diag,tmp)

    implicit none

    !> array with one batch of partial MO int.:
    real(realk), intent(inout) :: gmo(:)
    !> array containing the previous contributions
    !  to this MO int. batch, packed.
    type(tensor), intent(in) :: pack_gmo
    !> index corresponding to the current batch:
    integer, intent(in) :: tile
    !> dimensions of array:
    integer, intent(in) :: ntot, dimP, dimQ
    !> Diagonal block ?
    logical, intent(in) :: diag
    !> working array:
    real(realk), intent(inout) :: tmp(:)

    integer :: s, r, rs, sr, q, ibat1, ibat2, ipack, nnod
    integer(kind=long) :: ncopy

    ipack  = 1
    nnod = 1
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
#endif

    ! 1st case: current batch corresponds to diagonal block.
    if (diag) then

       ! get batch from pdm:
       ncopy = ntot*(ntot+1)*dimP*(dimP+1)/4

       if (nnod>1.and.pack_gmo%itype==TT_TILED_DIST) then
          call time_start_phase(PHASE_COMM)
          call tensor_get_tile(pack_gmo,tile,tmp,ncopy)
          call time_start_phase(PHASE_WORK)
       else if (nnod>1.and.pack_gmo%itype==TT_TILED_REPL) then
          call dcopy(ncopy,pack_gmo%ti(tile)%t,1,tmp,1)
       else
          call dcopy(ncopy,pack_gmo%elm2(1,tile),1,tmp,1)
       end if

       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ntot,gmo,tmp,dimP,dimQ)&
       !$OMP PRIVATE(r,s,rs,sr,ibat1,ibat2,ipack)
       do s=1,ntot
          do r=1,s
             rs = r + (s-1)*ntot
             do q=1,dimQ
                ipack = 1 + q*(q-1)/2 + (s*(s-1)/2 + r-1)*(dimQ*(dimQ+1)/2)
                ibat1 = 1 + (q-1)*dimP + (rs-1)*dimP*dimQ
                call dcopy(q,tmp(ipack),1,gmo(ibat1),1)
                ibat2 = q + (rs-1)*dimP*dimQ
                call dcopy(q-1,tmp(ipack),1,gmo(ibat2),dimP)
             end do
             if (r/=s) then
                sr = s + (r-1)*ntot
                sr = 1 + (sr-1)*dimP*dimQ
                rs = 1 + (rs-1)*dimP*dimQ
                call dcopy(dimP*dimQ,gmo(rs),1,gmo(sr),1)
             end if
          end do
       end do
       !$OMP END PARALLEL DO

       ! 2nd case: current batch corresponds to an upper diagonal block,
    else

       ! get batch from pdm:
       ncopy = dimP*dimQ*ntot*(ntot+1)/2

       if (nnod>1.and.pack_gmo%itype==TT_TILED_DIST) then
          call time_start_phase(PHASE_COMM)
          call tensor_get_tile(pack_gmo,tile,tmp,ncopy)
          call time_start_phase(PHASE_WORK)
       else if (nnod>1.and.pack_gmo%itype==TT_TILED_REPL) then
          call dcopy(ncopy,pack_gmo%ti(tile)%t,1,tmp,1)
       else
          call dcopy(ncopy,pack_gmo%elm2(1,tile),1,tmp,1)
       end if

       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ntot,gmo,tmp,dimP,dimQ)&
       !$OMP PRIVATE(r,s,rs,sr,ibat1,ipack,ncopy)
       do s=1,ntot
          do r=1,s
             rs = r + (s-1)*ntot
             ipack = 1 + (s*(s-1)/2 + r-1)*dimP*dimQ
             ibat1 = 1 + (rs-1)*dimP*dimQ
             ncopy = dimP*dimQ
             call dcopy(ncopy,tmp(ipack),1,gmo(ibat1),1)

             if (r/=s) then
                sr = s + (r-1)*ntot
                sr = 1 + (sr-1)*dimP*dimQ
                call dcopy(dimP*dimQ,gmo(ibat1),1,gmo(sr),1)
             end if
          end do
       end do
       !$OMP END PARALLEL DO

    end if

  end subroutine unpack_gmo
#endif

  subroutine get_mo_integral_par(integral,trafo1,trafo2,trafo3,trafo4,mylsitem,INTSPEC,local,collective)
    implicit none
    type(tensor),intent(inout)  :: integral
    type(tensor),intent(inout)  :: trafo1,trafo2,trafo3,trafo4
    type(lsitem), intent(inout) :: mylsitem
    character, intent(inout)    :: INTSPEC(5)
    logical, intent(in) :: local
    logical, intent(inout) :: collective
    !Integral stuff
    logical :: save_cs_screen, save_ps_screen
    integer :: alphaB,gammaB,dimAlpha,dimGamma
    integer :: dim1,dim2,dim3,MinAObatch
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
    integer :: iorb,nthreads
    integer :: idx,nb,n1,n2,n3,n4,fa,fg,la,lg,i,k,myload,nba,nbg,biA,biG,bsA,bsG
    logical :: FoundInMem,doscreen
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    Character(80)        :: FilenameCS,FilenamePS
    Character(80),pointer:: BatchfilenamesCS(:,:)
    Character(80),pointer:: BatchfilenamesPS(:,:)
    integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
    TYPE(DECscreenITEM)  :: DecScreen
    integer, pointer :: batchdimAlpha(:),batchdimGamma(:)
    logical :: fullRHS
    integer :: MaxAllowedDimAlpha,MaxActualDimAlpha,nbatchesAlpha
    integer :: MaxAllowedDimGamma,MaxActualDimGamma,nbatchesGamma
    real(realk), pointer :: w1(:),w2(:)
    real(realk) :: MemFree,nrm, MemToUse
    integer(kind=long) :: maxsize, addsize
    logical :: master
    integer(kind=ls_mpik) :: me, nnod
    integer               :: lenI2, modeBidx(2), modeBdim(2), el,jobidx,pos, bidx, bpos
    integer, pointer      :: jobdist(:)
    type(c_ptr)           :: jobdistc
    integer(kind=ls_mpik) :: jobdistw
    real(realk), pointer :: work(:)
    integer(kind=long) :: w1size, w2size
    real(realk), parameter :: fraction_of = 0.8E0_realk
    logical :: dynamic_load,completely_distributed, use_bg_buf, mem_saving
    integer :: nbuffs
    type(tensor) :: Cint, int1, int2, int3, t1_par, t2_par, t3_fa, t4_fg
    integer :: ndimA,  ndimB,  ndimC,  ndimD
    integer :: ndimAs, ndimBs, ndimCs, ndimDs
    integer :: startA, startB, startC, startD
    integer :: dims(4), tdim(4), starts(4), order4(4), buf_done, buf_sent, tilenr
    integer :: bs,as,gs,n1s,n2s,n3s,n4s, inc,b,e,m,n,t1,t2,t3,t4, scheme
    real(realk), pointer :: one(:),two(:),thr(:)
    integer              :: onen,  twon,  thrn
    integer(kind=8) :: nbu
    real(realk) :: tot_intloop,         tot_intloop_min,     tot_intloop_max
    real(realk) :: time_t4fg_tot_max,   time_t4fg_tot_min,   time_t4fg_tot,    time_t4fg 
    real(realk) :: time_t3fa_tot_max,   time_t3fa_tot_min,   time_t3fa_tot,    time_t3fa
    real(realk) :: time_int1_tot_max,   time_int1_tot_min,   time_int1_tot ,   time_int1
    real(realk) :: time_cont1_tot_max,  time_cont1_tot_min,  time_cont1_tot,   time_cont1 
    real(realk) :: time_cont2_tot_max,  time_cont2_tot_min,  time_cont2_tot,   time_cont2
    real(realk) :: time_cont3_tot_max,  time_cont3_tot_min,  time_cont3_tot,   time_cont3
    real(realk) :: time_cont4_tot_max,  time_cont4_tot_min,  time_cont4_tot,   time_cont4
    real(realk) :: phase_cntrs(nphases), total_integral_time
    real(realk) :: flushing_time, flushing_time_min, flushing_time_max 
    real(realk) :: unlock_time,   unlock_time_min,   unlock_time_max
    real(realk) :: waiting_time,  waiting_time_min,  waiting_time_max
    real(realk) :: time_w_min, time_w_max
    real(realk) :: time_c_min, time_c_max
    real(realk) :: time_i_min, time_i_max
    integer(kind=ls_mpik),pointer :: req(:)
#ifdef VAR_MPI
    integer(kind=ls_mpik), parameter :: mode = MPI_MODE_NOCHECK
#endif

    call time_start_phase( PHASE_WORK, twall = total_integral_time )

    master        = .true.
    me            = 0
    nnod          = 1
    dynamic_load  = DECinfo%dyn_load
    unlock_time   = 0.0E0_realk 
    waiting_time  = 0.0E0_realk
    flushing_time = 0.0E0_realk
#ifdef VAR_MPI
    master        = (infpar%lg_mynum == infpar%master)
    me            = infpar%lg_mynum
    nnod          = infpar%lg_nodtot
    unlock_time   = time_lsmpi_win_unlock 
    waiting_time  = time_lsmpi_wait
    flushing_time = time_lsmpi_win_flush
#endif
    completely_distributed = .false.
    if(.not.local)then
       nbuffs         = integral%ntpm(4)
    else
       nbuffs         = 2
    endif
    time_t4fg_tot  = 0.0E0_realk
    time_t3fa_tot  = 0.0E0_realk
    time_cont1_tot = 0.0E0_realk
    time_int1_tot  = 0.0E0_realk
    time_cont2_tot = 0.0E0_realk
    time_cont3_tot = 0.0E0_realk
    time_cont4_tot = 0.0E0_realk


    nb = trafo1%dims(1)
    n1 = trafo1%dims(2)
    n2 = trafo2%dims(2)
    n3 = trafo3%dims(2)
    n4 = trafo4%dims(2)
    !Get default splits, they should never be the ones used later
    n1s = get_split_scheme_0(n1)
    n2s = get_split_scheme_0(n2)
    n3s = get_split_scheme_0(n3)
    n4s = get_split_scheme_0(n4)

    !If the dim in the trafo matrix corresponds to the trafo dim then the
    !segmets are set to the same lengths. Same dimensions are required to
    !have the same segmenting, therefore it may be overwritten several times
    do i = 1, integral%mode
       if(integral%dims(i) == n1) n1s = integral%tdim(i)
       if(integral%dims(i) == n2) n2s = integral%tdim(i)
       if(integral%dims(i) == n3) n3s = integral%tdim(i)
       if(integral%dims(i) == n4) n4s = integral%tdim(i)
    enddo

    if( integral%dims(1) /= n1 .or. integral%dims(2) /= n2 .or. &
         & integral%dims(3) /= n3 .or. integral%dims(4) /= n4)then
       call lsquit("ERROR(get_mo_integral_par)wrong dimensions of the integrals&
            & or the transformation matrices",-1)
    endif
    bs = get_split_scheme_0(nb)

    IF(DECinfo%useIchor)THEN
       iprint     = 0       !print level for Ichor Integral code
       MoTrans    = .FALSE. !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol    = .TRUE.  !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       iAO = 1
       call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
       nullify(AOGammabatchinfo)
       nullify(AOalphabatchinfo)    
    ELSE
       doscreen                 = MyLsItem%setting%scheme%cs_screen.OR.MyLsItem%setting%scheme%ps_screen
       nullify(orb2batchAlpha)
       nullify(batchdimAlpha)
       nullify(batchsizeAlpha)
       nullify(batch2orbAlpha)
       nullify(batchindexAlpha)
       nullify(orb2batchGamma)
       nullify(batchdimGamma)
       nullify(batchsizeGamma)
       nullify(batch2orbGamma)
       nullify(batchindexGamma)
    ENDIF
#ifdef VAR_MPI
    if(master)then
       call time_start_phase( PHASE_COMM )
       if(.not.local)call wake_slaves_for_simple_mo(integral,trafo1,trafo2,trafo3,&
            &trafo4,mylsitem,collective)
       call time_start_phase( PHASE_WORK )
    endif
#endif
    !==================================================
    !                  Batch construction             !
    !==================================================

    use_bg_buf = mem_is_background_buf_init()
    if(use_bg_buf)then
       nbu = mem_get_bg_buf_free()
#ifdef VAR_MPI
       call lsmpi_reduce_min(nbu,infpar%master,infpar%lg_comm)
#endif
    else
       nbu = 0
    endif
    

    ! Get free memory and determine maximum batch sizes
    ! -------------------------------------------------
    if(master)then

       IF(DECinfo%useIchor)THEN
          !Determine the minimum allowed AObatch size MinAObatch
          !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
          !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
          !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
          !'R'  !Specifies that it is the Regular AO basis that should be batched
          iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
          call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
       ELSE
          call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
       ENDIF

       call get_currently_available_memory(MemFree)
       
       MemToUse = MemFree * fraction_of
       if( use_bg_buf )then
          MemToUse = MemToUse + (dble(nbu)*8.0E0_realk)/(1024.0E0_realk**3)
       endif

       if(DECinfo%PL>3)then
          write(*,'("CC integrals: operating with ",g9.2,"GB")')MemToUse
       endif

       call get_max_batch_and_scheme_ccintegral(maxsize,MinAObatch,scheme,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
       & n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,nbuffs,nbu,MyLsItem%setting,MemToUse,use_bg_buf)

       if(DECinfo%PL>2)then
          print *,"INFO(get_mo_integral_par): Requesting scheme:",scheme
          print *,"INFO(get_mo_integral_par): with Alpha:",MaxAllowedDimAlpha," Gamma:",MaxAllowedDimGamma
       endif
    endif


    if(.not.local)then
       integral%access_type = AT_ALL_ACCESS
       trafo1%access_type   = AT_ALL_ACCESS
       trafo2%access_type   = AT_ALL_ACCESS
       trafo3%access_type   = AT_ALL_ACCESS
       trafo4%access_type   = AT_ALL_ACCESS
#ifdef VAR_MPI
       call time_start_phase( PHASE_COMM )
       call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
       call ls_mpi_buffer(MaxAllowedDimAlpha,infpar%master)
       call ls_mpi_buffer(MaxAllowedDimGamma,infpar%master)
       call ls_mpi_buffer(scheme,infpar%master)
       call ls_mpi_buffer(nbuffs,infpar%master)
       call ls_mpi_buffer(INTSPEC,5,infpar%master)
       call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
       call time_start_phase( PHASE_WORK )
#endif
    endif


    select case(scheme)
    case(0)
       collective             = .true.
       mem_saving             = .false.
       completely_distributed = .false.
    case(1)
       collective             = .false.
       mem_saving             = .false.
       completely_distributed = .false.
    case(2)
       collective             = .false.
       mem_saving             = .true.
       completely_distributed = .false.
    case(3)
       collective             = .false.
       mem_saving             = .false.
       completely_distributed = .true.
    end select

    if(completely_distributed) then
       dynamic_load = .false.
       if( collective ) call lsquit("ERROR(get_mo_integral_par)&
       & completly_distributed and collective is impossible",-1)
    endif

    if(.not.completely_distributed)then
       ! ************************************************
       ! * Determine batch information for Gamma batch  *
       ! ************************************************
       IF(DECinfo%useIchor)THEN
          iAO = 4 !Gamma is the 4. Center of the 4 center two electron coulomb integral
          !Determine how many batches of AOS based on the MaxAllowedDimGamma, the requested
          !size of the AO batches. iAO is the center that the batching should occur on. 
          !'R'  !Specifies that it is the Regular AO basis that should be batched 
          call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
               & nbatchesGamma,DECinfo%output)
          call mem_alloc(AOGammabatchinfo,nbatchesGamma)
          !Construct the batches of AOS based on the MaxAllowedDimGamma, the requested
          !size of the AO batches - MaxAllowedDimGamma must be unchanged since the call 
          !to determine_Ichor_nbatchesofAOS
          !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
          !So MaxActualDimGamma must be less og equal to MaxAllowedDimGamma
          call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
               & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
       ELSE
          ! Orbital to batch information
          ! ----------------------------
          call mem_alloc(orb2batchGamma,nb)
          call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
               & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
               &nbatchesGamma,orb2BatchGamma,'R')
       ENDIF

       if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma,&
          & 'with maximum size',MaxActualDimGamma

       IF(.NOT.DECinfo%useIchor)THEN
          ! Translate batchindex to orbital index
          ! -------------------------------------
          call mem_alloc(batch2orbGamma,nbatchesGamma)
          do idx=1,nbatchesGamma
             call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx))
             batch2orbGamma(idx)%orbindex = 0
             batch2orbGamma(idx)%norbindex = 0
          end do
          do iorb=1,nb
             idx = orb2batchGamma(iorb)
             batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
             K = batch2orbGamma(idx)%norbindex
             batch2orbGamma(idx)%orbindex(K) = iorb
          end do
       ENDIF

       ! ************************************************
       ! * Determine batch information for Alpha batch  *
       ! ************************************************

       IF(DECinfo%useIchor)THEN
          iAO = 3 !Alpha is the 3. Center of the 4 center two electron coulomb integral
          !Determine how many batches of AOS based on the MaxAllowedDimAlpha, the requested
          !size of the AO batches. iAO is the center that the batching should occur on. 
          !'R'  !Specifies that it is the Regular AO basis that should be batched 
          call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
               & nbatchesAlpha,DECinfo%output)
          call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
          !Construct the batches of AOS based on the MaxAllowedDimAlpha, the requested
          !size of the AO batches - MaxAllowedDimAlpha must be unchanged since the call 
          !to determine_Ichor_nbatchesofAOS
          !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
          !So MaxActualDimAlpha must be less og equal to MaxAllowedDimAlpha
          call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
               & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
       ELSE
          ! Orbital to batch information
          ! ----------------------------
          call mem_alloc(orb2batchAlpha,nb)
          call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
               & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')
       ENDIF

       if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha&
          &, 'with maximum size',MaxActualDimAlpha

       IF(.NOT.DECinfo%useIchor)THEN
          ! Translate batchindex to orbital index
          ! -------------------------------------
          call mem_alloc(batch2orbAlpha,nbatchesAlpha)
          do idx=1,nbatchesAlpha
             call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
             batch2orbAlpha(idx)%orbindex = 0
             batch2orbAlpha(idx)%norbindex = 0
          end do
          do iorb=1,nb
             idx = orb2batchAlpha(iorb)
             batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
             K = batch2orbAlpha(idx)%norbindex
             batch2orbAlpha(idx)%orbindex(K) = iorb
          end do
       ENDIF
    else

       MaxActualDimAlpha = MaxAllowedDimAlpha
       MaxActualDimGamma = MaxAllowedDimGamma
       save_cs_screen = mylsitem%setting%SCHEME%CS_SCREEN
       save_ps_screen = mylsitem%setting%SCHEME%PS_SCREEN
       mylsitem%setting%SCHEME%CS_SCREEN = .FALSE.
       mylsitem%setting%SCHEME%PS_SCREEN = .FALSE.
       doscreen = mylsitem%setting%SCHEME%CS_SCREEN.OR.mylsitem%setting%SCHEME%PS_SCREEN


       call tensor_ainit( t1_par, [nb,n1], 2, local=local, atype="TDAR", tdims=[bs,n1s] )
       call tensor_convert(trafo1%elm1,t1_par)
       call tensor_ainit( t2_par, [nb,n2], 2, local=local, atype="TDAR", tdims=[bs,n2s] )
       call tensor_convert(trafo2%elm1,t2_par)

       nbatchesAlpha = nb / MaxActualDimAlpha
       if( mod( nb, MaxActualDimAlpha ) > 0 ) nbatchesAlpha = nbatchesAlpha + 1
       nbatchesGamma = nb / MaxActualDimGamma
       if( mod( nb, MaxActualDimGamma ) > 0 ) nbatchesGamma = nbatchesGamma + 1

    endif

    w1size=get_work_array_size(1,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,MaxActualDimAlpha,MaxActualDimGamma,&
       & scheme,nbuffs,MyLsItem%setting)
    w2size=get_work_array_size(2,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,MaxActualDimAlpha,MaxActualDimGamma,&
       & scheme,nbuffs,MyLsItem%setting)

    onen = (n1s*nb)*MaxActualDimAlpha*MaxActualDimGamma
    twon = (n1s*n2s)*MaxActualDimAlpha*MaxActualDimGamma
    thrn = (n1s*n2s)*n3s*MaxActualDimGamma

    maxsize = w1size + w2size
    addsize = 0
    if( collective ) addsize = (i8*n1)*n2*n3*n4
    if( mem_saving ) addsize = nbuffs*(i8*n1s)*n2s*n3s*n4s
    maxsize = maxsize + addsize
    
    if(master)then
       print *,"INFO(get_mo_integral_par): Getting Alpha:",MaxActualDimGamma," Gamma:",MaxActualDimGamma
       print *,"INFO(get_mo_integral_par): with elements:",maxsize,"in buffer",nbu
    endif

    if( use_bg_buf ) then
       if(maxsize > nbu) then
          print *, "Warning(get_mo_integral_par):  This should not happen, if the memory counting is correct&
             &, Node:",me," requests ",maxsize," in buffer ",nbu
          !call mem_change_background_alloc(maxsize*8_long)
       endif

       call mem_pseudo_alloc( w1, w1size )
       call mem_pseudo_alloc( w2, w2size )
    else
       call mem_alloc( w1, w1size )
       call mem_alloc( w2, w2size )
    endif

    !First touch
    w1(1) = 0.0E0_realk
    w2(1) = 0.0E0_realk

    if(collective.or.mem_saving)then

       if(mem_saving)call mem_alloc(req,nbuffs)

       if( use_bg_buf ) then
          call mem_pseudo_alloc(work,addsize)
       else
          call mem_alloc(work,addsize)
       endif

       work = 0.0E0_realk

    endif


    ! ************************************************
    ! *  precalculate the full schreening matrix     *
    ! ************************************************
    if(.not.completely_distributed)then
       IF(DECinfo%useIchor)THEN
          !Calculate Screening integrals 
          SameMOL = .TRUE. !Specifies same molecule on all centers 
          call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
       ELSE
          ! This subroutine builds the full screening matrix.
          call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
               & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
          if(mylsitem%setting%scheme%cs_screen .OR. mylsitem%setting%scheme%ps_screen)THEN
             call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
                  & nb,nbatchesAlpha,nbatchesGamma,&
                  & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
                  & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
             call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
                  & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
                  & batchindexAlpha,batchindexGamma,&
                  & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
          endif
       ENDIF
    endif


#ifdef VAR_MPI
    if(.not.completely_distributed)then
       IF(DECinfo%useIchor)THEN
          call mem_alloc(batchdimAlpha,nbatchesAlpha)
          do idx=1,nbatchesAlpha
             batchdimAlpha(idx) = AOAlphabatchinfo(idx)%dim 
          enddo
          call mem_alloc(batchdimGamma,nbatchesGamma)
          do idx=1,nbatchesGamma
             batchdimGamma(idx) = AOGammabatchinfo(idx)%dim 
          enddo
       ENDIF
       !JOB distribution
       if(.not.dynamic_load)then
          lenI2 = nbatchesAlpha*nbatchesGamma
          call mem_alloc(jobdist,lenI2)
          myload   = 0
          jobdist  = 0
          call distribute_mpi_jobs(jobdist,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
             &batchdimGamma,myload,nnod,me)
       else
          lenI2 = 1
          call mem_alloc( jobdist, jobdistc, lenI2 ) 
          jobdist = 0
          if(me == infpar%master) jobdist(1) = nnod + 1

          call lsmpi_win_create(jobdist,jobdistw,1,infpar%lg_comm)
#ifdef VAR_HAVE_MPI3
          call lsmpi_win_lock_all(jobdistw,ass=mode)
#endif

       endif
       IF(DECinfo%useIchor)THEN
          call mem_dealloc(batchdimAlpha)
          call mem_dealloc(batchdimGamma)
       ENDIF
    else
       !ALL nodes go through the loops together
       lenI2 = nbatchesAlpha*nbatchesGamma
       call mem_alloc(jobdist,lenI2)
       myload   = 0
       jobdist  = me
    endif
#else
    call mem_alloc(jobdist,nbatchesAlpha*nbatchesGamma)
    jobdist = 0
#endif

    myload = 0
    fullRHS = nbatchesGamma.EQ.1.AND.nbatchesAlpha.EQ.1

#ifdef VAR_MPI
    if(.not.collective .and. alloc_in_dummy)then
       call tensor_lock_wins(integral, 's', all_nodes = .true. )
    endif
#endif


    gammaB   = 0
    alphaB   = 0
    modeBdim = [nbatchesAlpha,nbatchesGamma]
    jobidx   = 0
    buf_sent = 0
    buf_done = 0

    call time_start_phase( PHASE_WORK, twall = tot_intloop )
    call time_phases_get_current(current_wt=phase_cntrs)

    BatchLoop: do while(gammaB <= nbatchesGamma.or.alphaB<= nbatchesAlpha)  ! AO batches


       !SET ALPHAB AND GAMMAB
       if(.not.dynamic_load)then

          jobidx = jobidx + 1

          if( jobidx > nbatchesGamma * nbatchesAlpha ) exit BatchLoop

          call get_midx(jobidx,modeBidx,modeBdim,2)

          gammaB = modeBidx(2)
          alphaB = modeBidx(1)

          if( me /= jobdist(gammaB + (alphaB-1) *nbatchesGamma) ) cycle BatchLoop

       else

          if(alphaB == 0)then
             call get_midx(int(me + 1),modeBidx,modeBdim,2)
             gammaB = modeBidx(2)
             alphaB = modeBidx(1)
          else
#ifdef VAR_MPI
             el  = 1
             pos = 1
#ifdef VAR_HAVE_MPI3
             call lsmpi_get_acc(el,jobidx,infpar%master,pos,jobdistw)
             call lsmpi_win_flush(jobdistw,rank = infpar%master, local=.true.)
#else
             call lsmpi_win_lock(infpar%master,jobdistw,'e')
             call lsmpi_get_acc(el,jobidx,infpar%master,pos,jobdistw)
             call lsmpi_win_unlock(infpar%master,jobdistw)
#endif
#endif

             if( jobidx > nbatchesGamma * nbatchesAlpha ) exit BatchLoop

             call get_midx(jobidx,modeBidx,modeBdim,2)
             gammaB = modeBidx(2)
             alphaB = modeBidx(1)

          endif



       endif


       if(alphaB > nbatchesAlpha .or. gammaB > nbatchesGamma)then
          write (*, '("Rank",I3," has invalid job (",I3,"/",I3,",",I3,"/",I3,")")') &
             & me,alphaB,nbatchesAlpha,gammaB,nbatchesGamma
          call lsquit("ERROR(get_mo_integral_par) error in jobassignment",-1)
       endif


       if(DECinfo%PL>2)then
          if(completely_distributed)then
             if( me == 0 )then
                write (*, '("Starting job (",I3,"/",I3,",",I3,"/",I3,")")') &
                   & alphaB,nbatchesAlpha,gammaB,nbatchesGamma
             endif
          else
             write (*, '("Rank",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")') &
                & me,alphaB,nbatchesAlpha,gammaB,nbatchesGamma
          endif
       endif

       if(DECinfo%ccsolverskip)then
          call random_number(w1)
          call random_number(w2)
       else
       if(.not.completely_distributed)then
          IF(DECinfo%useIchor)THEN
             lg = AOGammabatchinfo(gammaB)%dim               ! Dimension of gamma batch
             fg = AOGammabatchinfo(gammaB)%orbstart          ! First orbital index in gamma batch
             GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
             AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
             AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
          ELSE
             lg  = batchdimGamma(gammaB)                     ! Dimension of gamma batch
             fg  = batch2orbGamma(gammaB)%orbindex(1)        ! First index in gamma batch
             biG = batchindexGamma(gammaB)
             bsG = batchsizeGamma(gammaB)
          ENDIF

          IF(DECinfo%useIchor)THEN
             la = AOAlphabatchinfo(alphaB)%dim               ! Dimension of alpha batch
             fa = AOAlphabatchinfo(alphaB)%orbstart          ! First orbital index in alpha batch
             AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
             AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
             AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
          ELSE
             la  = batchdimAlpha(alphaB)                              ! Dimension of alpha batch
             fa  = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
             biA = batchindexAlpha(alphaB)
             bsA = batchsizeAlpha(alphaB)
          ENDIF
          !print '(I3,"have",8I7)',me,lg,fg,biG,bsG,la,fa,biA,bsA
          !call lsmpi_barrier(infpar%lg_comm)

          myload     = myload + la * lg


          call time_start_phase(PHASE_WORK, twall = time_int1 )
          IF(DECinfo%useIchor)THEN
             call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nb,nb,la,lg,&
                  & w1,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
                  & AOGammaStart,AOGammaEnd,MoTrans,nb,nb,la,lg,NoSymmetry,DECinfo%IntegralThreshold)
          ELSE
             IF(doscreen) Mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
             IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p
             
             call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, Mylsitem%setting, w1,biA,&
                  &biG,bsA,bsG,nb,nb,la,lg,fullRHS,INTSPEC,DECinfo%IntegralThreshold)
          ENDIF
          call time_start_phase(PHASE_WORK, ttot = time_int1 )
          time_int1_tot = time_int1_tot + time_int1


#ifdef VAR_MPI
          if( .not.(collective.or.mem_saving).and.alloc_in_dummy )then
             call lsmpi_win_flush(integral%wi(1),local=.true.)
          endif
#endif



          if(mem_saving)then

#ifdef VAR_MPI


             one => w2(1:onen)
             two => w2(onen+1:onen+twon)
             thr => w2(onen+twon+1:onen+twon+thrn)

             !TODO: introduce OMP in order to speed up the routine
             do t1 = 1, integral%ntpm(1)

                startA = 1 + (t1-1)*integral%tdim(1)
                if(((n1-(t1-1)*n1s)/n1s)>=1)then
                   ndimA=n1s
                else
                   ndimA=mod(n1,n1s)
                endif

                m = nb*la*lg
                n = ndimA
                k = nb
                call time_start_phase(PHASE_WORK, twall = time_cont1 )
                call dgemm('t','n',m,n,k,1.0E0_realk,w1,k,trafo1%elm2(1,startA),nb,0.0E0_realk,one,m)
                call time_start_phase(PHASE_WORK, ttot = time_cont1 )
                time_cont1_tot = time_cont1_tot + time_cont1

                do t2 = 1, integral%ntpm(2)

                   startB = 1 + (t2-1)*integral%tdim(2)
                   if(((n2-(t2-1)*n2s)/n2s)>=1)then
                      ndimB=n2s
                   else
                      ndimB=mod(n2,n2s)
                   endif

                   m = la*lg*ndimA
                   n = ndimB
                   k = nb
                   call time_start_phase(PHASE_WORK, twall = time_cont2 )
                   call dgemm('t','n',m,n,k,1.0E0_realk,one,k,trafo2%elm2(1,startB),nb,0.0E0_realk,two,m)
                   call time_start_phase(PHASE_WORK, ttot = time_cont2 )
                   time_cont2_tot = time_cont2_tot + time_cont2

                   do t3 = 1, integral%ntpm(3)

                      startC = 1 + (t3-1)*integral%tdim(3)
                      if(((n3-(t3-1)*n3s)/n3s)>=1)then
                         ndimC=n3s
                      else
                         ndimC=mod(n3,n3s)
                      endif
                      m = lg*ndimA*ndimB
                      n = ndimC
                      k = la

                      call time_start_phase(PHASE_WORK, twall = time_cont3 )
                      call dgemm('t','n',m,n,k,1.0E0_realk,two,k,trafo3%elm2(fa,startC),nb,0.0E0_realk,thr,m)
                      call time_start_phase(PHASE_WORK, ttot = time_cont3 )
                      time_cont3_tot = time_cont3_tot + time_cont3

                      do t4 = 1, integral%ntpm(4)

                         startD = 1 + (t4-1)*integral%tdim(4)
                         if(((n4-(t4-1)*n4s)/n4s)>=1)then
                            ndimD=n4s
                         else
                            ndimD=mod(n4,n4s)
                         endif

                         m = ndimA*ndimB*ndimC
                         n = ndimD
                         k = lg


                         tilenr=get_cidx([t1,t2,t3,t4],integral%ntpm,integral%mode)

                         bidx = mod(buf_sent,nbuffs)+1

                         bpos = 1 + (bidx-1) * n1s*n2s*n3s*n4s

                         if(alloc_in_dummy .and. buf_sent>=nbuffs)then
                            call lsmpi_wait(req(bidx))
                            buf_done = buf_done + 1
                         endif
    

                         call dgemm('t','n',m,n,k,1.0E0_realk,thr,k,trafo4%elm2(fg,startD),nb,0.0E0_realk,work(bpos),m)
                         !accumulate TODO, meeds to be optimized
#ifdef VAR_HAVE_MPI3
                         call tensor_accumulate_tile(integral,[t1,t2,t3,t4],work(bpos:bpos+m*n-1),m*n,&
                            &lock_set=.true.,req=req(bidx))
#else
                         call tensor_accumulate_tile(integral,tilenr,work(bpos:bpos+m*n-1),m*n,&
                            &lock_set=.false.)
#endif
                         buf_sent = buf_sent + 1
                      enddo
                   enddo
                enddo
             enddo
#else
             call lsquit("ERROR(get_mo_integral_par):this option is MPI only",-1)
#endif


          else
             ndimA  = n1
             ndimB  = n2
             ndimC  = n3
             ndimD  = n4

             startA = 1
             startB = 1
             startC = 1
             startD = 1

             one => w2
             two => w1
             thr => w2

             !something more sophisticated can be implemented here
             m = nb*la*lg
             n = ndimA
             k = nb
             call time_start_phase(PHASE_WORK, twall = time_cont1 )
             call dgemm('t','n',m,n,k,1.0E0_realk,w1,k,trafo1%elm2(1,startA),nb,0.0E0_realk,one,m)
             call time_start_phase(PHASE_WORK, ttot = time_cont1 )
             time_cont1_tot = time_cont1_tot + time_cont1

             m = la*lg*ndimA
             n = ndimB
             k = nb
             call time_start_phase(PHASE_WORK, twall = time_cont2 )
             call dgemm('t','n',m,n,k,1.0E0_realk,one,k,trafo2%elm2(1,startB),nb,0.0E0_realk,two,m)
             call time_start_phase(PHASE_WORK, ttot = time_cont2 )
             time_cont2_tot = time_cont2_tot + time_cont2

             m = lg*ndimA*ndimB
             n = ndimC
             k = la
             call time_start_phase(PHASE_WORK, twall = time_cont3 )
             call dgemm('t','n',m,n,k,1.0E0_realk,two,k,trafo3%elm2(fa,startC),nb,0.0E0_realk,thr,m)
             call time_start_phase(PHASE_WORK, ttot = time_cont3 )
             time_cont3_tot = time_cont3_tot + time_cont3

             m = ndimA*ndimB*ndimC
             n = ndimD
             k = lg

             call time_start_phase(PHASE_WORK, twall = time_cont4 )
             if(collective) then
                call dgemm('t','n',m,n,k,1.0E0_realk,thr,k,trafo4%elm1(fg),nb,1.0E0_realk,work,m)
             else
                call dgemm('t','n',m,n,k,1.0E0_realk,thr,k,trafo4%elm1(fg),nb,0.0E0_realk,w1,m)

                call time_start_phase( PHASE_COMM )
                call tensor_add(integral,1.0E0_realk,w1,wrk=w2,iwrk=w2size)
                call time_start_phase( PHASE_WORK )
             endif
             call time_start_phase(PHASE_WORK, ttot = time_cont4 )
             time_cont4_tot = time_cont4_tot + time_cont4
          endif



       else

          !short hand notation
          fg = 1 + (gammaB-1)*MaxActualDimGamma
          lg = nb - fg + 1
          if( lg >= MaxActualDimGamma )then
             lg = MaxActualDimGamma
          endif

          if( lg > bs )then
             gs = bs
          else
             gs = lg
          endif
          
          call time_start_phase(PHASE_WORK, twall = time_t4fg )
          call tensor_ainit(t4_fg, [lg,n4], 2, local=local, atype="TDAR", tdims=[gs,n4s])
          call copy_stripe_from_full_matrix(trafo4%elm1,w1,fg,lg,nb,n4)
          call tensor_convert(w1,t4_fg)
          call time_start_phase(PHASE_WORK, ttot = time_t4fg )
          time_t4fg_tot = time_t4fg_tot + time_t4fg

          !short hand notation
          fa = 1 + (alphaB-1)*MaxActualDimAlpha
          la = nb - fa + 1
          if( la >= MaxActualDimAlpha )then
             la = MaxActualDimAlpha
          endif

          if( la > bs )then
             as = bs
          else
             as = la
          endif

          myload = myload + la * lg

          call time_start_phase(PHASE_WORK, twall = time_t3fa )
          call tensor_ainit(t3_fa, [la,n3], 2, local=local, atype="TDAR", tdims=[as,n3s])
          call copy_stripe_from_full_matrix(trafo3%elm1,w1,fa,la,nb,n3)
          call tensor_convert(w1,t3_fa)
          call time_start_phase(PHASE_WORK, ttot = time_t3fa )
          time_t3fa_tot = time_t3fa_tot + time_t3fa

          call time_start_phase(PHASE_WORK, twall = time_int1 )
          call tensor_ainit(Cint, [nb,nb,la,lg], 4, local=local, atype="TDAR", tdims=[bs,bs,as,gs])
          do i = 1, Cint%nlti

             call get_midx(Cint%ti(i)%gt,starts,Cint%ntpm,Cint%mode)

             ndimA  = Cint%ti(i)%d(1)
             ndimB  = Cint%ti(i)%d(2)
             ndimC  = Cint%ti(i)%d(3)
             ndimD  = Cint%ti(i)%d(4)

             startA = 1  + (starts(1)-1)*Cint%tdim(1)
             startB = 1  + (starts(2)-1)*Cint%tdim(2)
             startC = fa + (starts(3)-1)*Cint%tdim(3)
             startD = fg + (starts(4)-1)*Cint%tdim(4)

             call II_GET_ERI_INTEGRALBLOCK_INQUIRE(DECinfo%output,DECinfo%output,Mylsitem%setting,&
                & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
                & ndimAs,ndimBs,ndimCs,ndimDs,INTSPEC)

             call II_GET_ERI_INTEGRALBLOCK(DECinfo%output,DECinfo%output,Mylsitem%setting,&
                  & startA,startB,startC,startD,ndimA,ndimB,ndimC,ndimD,&
                  & ndimAs,ndimBs,ndimCs,ndimDs,INTSPEC,Cint%ti(i)%t,w1,&
                  & DECinfo%IntegralThreshold,DECinfo%useIchor)

          enddo
          call time_start_phase(PHASE_WORK, ttot = time_int1 )
          time_int1_tot = time_int1_tot + time_int1

          !call print_norm(Cint,"Integral",print_on_rank=0)

          call time_start_phase(PHASE_WORK, twall = time_cont1 )
          call tensor_ainit( int1, [nb,la,lg,n1], 4, local=local, atype="TDAR", tdims=[bs,as,gs,n1s] )
          order4 = [2,3,4,1]
          call tensor_contract( 1.0E0_realk, t1_par,Cint,[1],[1],1,0.0E0_realk, int1, order4, &
             & force_sync=.true.,wrk=w1,iwrk=w1size)
          call time_start_phase(PHASE_WORK, ttot = time_cont1 )
          time_cont1_tot = time_cont1_tot + time_cont1

          !call print_norm(int1,"int1",print_on_rank=0)

          call time_start_phase(PHASE_WORK, twall = time_cont2 )
          call tensor_ainit( int2, [la,lg,n1,n2], 4, local=local, atype="TDAR", tdims=[as,gs,n1s,n2s] )
          call tensor_free( Cint )
          order4 = [2,3,4,1]
          call tensor_contract( 1.0E0_realk, t2_par,int1,[1],[1],1,0.0E0_realk, int2, order4, &
             & force_sync=.true.,wrk=w1,iwrk=w1size)
          call time_start_phase(PHASE_WORK, ttot = time_cont2 )
          time_cont2_tot = time_cont2_tot + time_cont2

          !call print_norm(int2,"int2",print_on_rank=0)

          call time_start_phase(PHASE_WORK, twall = time_cont3 )
          call tensor_ainit( int3, [lg,n1,n2,n3], 4, local=local, atype="TDAR", tdims=[gs,n1s,n2s,n3s] )
          call tensor_free(int1)
          order4 = [2,3,4,1]
          call tensor_contract( 1.0E0_realk, t3_fa,int2,[1],[1],1,0.0E0_realk, int3, order4, &
             & force_sync=.true.,wrk=w1,iwrk=w1size)
          call time_start_phase(PHASE_WORK, ttot = time_cont3 )
          time_cont3_tot = time_cont3_tot + time_cont3

          !call print_norm(int3,"int3",print_on_rank=0)

          call time_start_phase(PHASE_WORK, twall = time_cont4 )
          call tensor_free(int2)
          order4 = [2,3,4,1]
          call tensor_contract( 1.0E0_realk, t4_fg,int3,[1],[1],1,1.0E0_realk, integral, order4, &
             & force_sync=.true.,wrk=w1,iwrk=w1size)
          call time_start_phase(PHASE_WORK, ttot = time_cont4 )
          time_cont4_tot = time_cont4_tot + time_cont4

          !call print_norm(integral,"integral",print_on_rank=0)

          call tensor_free( int3 )

          call tensor_free( t3_fa )
          call tensor_free( t4_fg )
       endif
       endif

    enddo BatchLoop

    call time_phases_get_diff(current_wt=phase_cntrs)
    call time_start_phase( PHASE_WORK, ttot = tot_intloop )

    ! Free integral stuff
    ! *******************
    if( .not. completely_distributed )then

       IF(DECinfo%useIchor)THEN
          call FREE_SCREEN_ICHORERI()
          call mem_dealloc(AOGammabatchinfo)
          call mem_dealloc(AOAlphabatchinfo)
       ELSE
          nullify(Mylsitem%setting%LST_GAB_LHS)
          nullify(Mylsitem%setting%LST_GAB_RHS)
          call free_decscreen(DECSCREEN)
          
          ! Free gamma stuff
          call mem_dealloc(orb2batchGamma)
          call mem_dealloc(batchdimGamma)
          call mem_dealloc(batchsizeGamma)
          call mem_dealloc(batchindexGamma)
          do i=1,nbatchesGamma
             call mem_dealloc(batch2orbGamma(i)%orbindex)
             batch2orbGamma(i)%orbindex => null()
          end do
          call mem_dealloc(batch2orbGamma)
          
          ! Free alpha stuff
          call mem_dealloc(orb2batchAlpha)
          call mem_dealloc(batchdimAlpha)
          call mem_dealloc(batchsizeAlpha)
          call mem_dealloc(batchindexAlpha)
          do i=1,nbatchesAlpha
             call mem_dealloc(batch2orbAlpha(i)%orbindex)
             batch2orbAlpha(i)%orbindex => null()
          end do
          call mem_dealloc(batch2orbAlpha)
       ENDIF

    else

       mylsitem%setting%SCHEME%CS_SCREEN = save_cs_screen
       mylsitem%setting%SCHEME%PS_SCREEN = save_ps_screen
       call tensor_free( t1_par )
       call tensor_free( t2_par )

    endif


#ifdef VAR_MPI
    if(.not.collective .and. alloc_in_dummy)then
       if(mem_saving)then

          do t1 = buf_done+1, buf_sent
             bidx = mod(t1-1,nbuffs)+1
             call lsmpi_wait(req(bidx))
          enddo

          call mem_dealloc(req)

       endif

       call tensor_unlock_wins(integral, all_nodes = .true. )

    endif

    call time_start_phase( PHASE_IDLE )
    call lsmpi_barrier(infpar%lg_comm)

    if(collective)then
       call time_start_phase( PHASE_COMM )
       call lsmpi_allreduce(work,(i8*n1)*n2*n3*n4,infpar%lg_comm)
       call time_start_phase( PHASE_WORK )
       call tensor_convert(work,integral, wrk = w1, iwrk = w1size   )
    endif

#else

    call time_start_phase( PHASE_WORK )
    call tensor_convert(work,integral, wrk = w1, iwrk = w1size )

#endif


    if(collective.or.mem_saving)then
       if( use_bg_buf )then
          call mem_pseudo_dealloc( work )
       else
          call mem_dealloc( work )
       endif
    endif
    if( use_bg_buf )then
       call mem_pseudo_dealloc( w2 )
       call mem_pseudo_dealloc( w1 )
    else
       call mem_dealloc( w2 )
       call mem_dealloc( w1 )
    endif

    if(DECinfo%ccsolverskip)then
       call tensor_random(integral)
    endif


    if(.not.dynamic_load)then
       call mem_dealloc(jobdist)
    else
#ifdef VAR_MPI
#ifdef VAR_HAVE_MPI3
       call lsmpi_win_unlock_all(jobdistw)
#endif
       call lsmpi_win_free(jobdistw)
       call mem_dealloc(jobdist,jobdistc)
#endif
    endif

    !TIMING INFORMATION
    if(DECinfo%PL>2)then

       tot_intloop_min    = tot_intloop 
       time_t4fg_tot_min  = time_t4fg_tot   
       time_t3fa_tot_min  = time_t3fa_tot   
       time_int1_tot_min  = time_int1_tot   
       time_cont1_tot_min = time_cont1_tot   
       time_cont2_tot_min = time_cont2_tot   
       time_cont3_tot_min = time_cont3_tot   
       time_cont4_tot_min = time_cont4_tot   

#ifdef VAR_MPI
       call lsmpi_reduce_realk_min( tot_intloop_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_t4fg_tot_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_t3fa_tot_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_int1_tot_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_cont1_tot_min, infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_cont2_tot_min, infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_cont3_tot_min, infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_cont4_tot_min, infpar%master, infpar%lg_comm )
#endif

       tot_intloop_max    = tot_intloop 
       time_t4fg_tot_max  = time_t4fg_tot   
       time_t3fa_tot_max  = time_t3fa_tot   
       time_int1_tot_max  = time_int1_tot   
       time_cont1_tot_max = time_cont1_tot   
       time_cont2_tot_max = time_cont2_tot   
       time_cont3_tot_max = time_cont3_tot   
       time_cont4_tot_max = time_cont4_tot   

#ifdef VAR_MPI
       call lsmpi_reduce_realk_max( tot_intloop_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_t4fg_tot_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_t3fa_tot_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_int1_tot_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_cont1_tot_max, infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_cont2_tot_max, infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_cont3_tot_max, infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_cont4_tot_max, infpar%master, infpar%lg_comm )

       call lsmpi_local_reduction( tot_intloop   , infpar%master )
       call lsmpi_local_reduction( time_t4fg_tot , infpar%master )
       call lsmpi_local_reduction( time_t3fa_tot , infpar%master )
       call lsmpi_local_reduction( time_int1_tot , infpar%master )
       call lsmpi_local_reduction( time_cont1_tot, infpar%master )
       call lsmpi_local_reduction( time_cont2_tot, infpar%master )
       call lsmpi_local_reduction( time_cont3_tot, infpar%master )
       call lsmpi_local_reduction( time_cont4_tot, infpar%master )

       call lsmpi_local_reduction(phase_cntrs,nphases,infpar%master)
#endif
       unlock_time   = time_lsmpi_win_unlock - unlock_time
       waiting_time  = time_lsmpi_wait       - waiting_time
       flushing_time = time_lsmpi_win_flush  - flushing_time

       unlock_time_min    = unlock_time
       waiting_time_min   = waiting_time      
       flushing_time_min  = flushing_time 

       unlock_time_max    = unlock_time
       waiting_time_max   = waiting_time      
       flushing_time_max  = flushing_time 

#ifdef VAR_MPI
       call lsmpi_reduce_realk_min( unlock_time_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( waiting_time_min  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( flushing_time_min , infpar%master, infpar%lg_comm )

       call lsmpi_reduce_realk_max( unlock_time_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( waiting_time_max  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( flushing_time_max , infpar%master, infpar%lg_comm )

       call lsmpi_local_reduction( unlock_time   , infpar%master )
       call lsmpi_local_reduction( waiting_time  , infpar%master )
       call lsmpi_local_reduction( flushing_time , infpar%master )
#endif

       time_w_min = phase_cntrs( PHASE_WORK_IDX )
       time_c_min = phase_cntrs( PHASE_COMM_IDX )
       time_i_min = phase_cntrs( PHASE_IDLE_IDX )

       time_w_max = phase_cntrs( PHASE_WORK_IDX )
       time_c_max = phase_cntrs( PHASE_COMM_IDX )
       time_i_max = phase_cntrs( PHASE_IDLE_IDX )

#ifdef VAR_MPI
       call lsmpi_reduce_realk_min( time_w_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_c_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_i_min , infpar%master, infpar%lg_comm )

       call lsmpi_reduce_realk_max( time_w_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_c_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_i_max , infpar%master, infpar%lg_comm )

       call lsmpi_local_reduction(phase_cntrs,nphases,infpar%master)
#endif

       if(master)then
          write(*,'("INTEGRAL total time               ",g10.3,g10.3,g10.3)')&
             & tot_intloop_max     ,      tot_intloop      /dble(nnod),tot_intloop_min 
          write(*,'("INTEGRAL time_t4fg_tot           ",g10.3,g10.3,g10.3,g10.3)')&
             & time_t4fg_tot_max  ,  time_t4fg_tot   /dble(nnod),time_t4fg_tot_min  ,  time_t4fg_tot   / tot_intloop
          write(*,'("INTEGRAL time_t3fa_tot           ",g10.3,g10.3,g10.3,g10.3)')&
             & time_t3fa_tot_max  ,  time_t3fa_tot   /dble(nnod),time_t3fa_tot_min  ,  time_t3fa_tot   / tot_intloop
          write(*,'("INTEGRAL time_int1_tot            ",g10.3,g10.3,g10.3,g10.3)')&
             & time_int1_tot_max   ,  time_int1_tot    /dble(nnod),time_int1_tot_min   ,  time_int1_tot    / tot_intloop
          write(*,'("INTEGRAL time_cont1_tot           ",g10.3,g10.3,g10.3,g10.3)')&
             & time_cont1_tot_max  ,  time_cont1_tot   /dble(nnod),time_cont1_tot_min  ,  time_cont1_tot   / tot_intloop
          write(*,'("INTEGRAL time_cont2_tot           ",g10.3,g10.3,g10.3,g10.3)')&
             & time_cont2_tot_max  ,  time_cont2_tot   /dble(nnod),time_cont2_tot_min  ,  time_cont2_tot   / tot_intloop
          write(*,'("INTEGRAL time_cont3_tot           ",g10.3,g10.3,g10.3,g10.3)')&
             & time_cont3_tot_max  ,  time_cont3_tot   /dble(nnod),time_cont3_tot_min  ,  time_cont3_tot   / tot_intloop
          write(*,'("INTEGRAL time_cont4_tot           ",g10.3,g10.3,g10.3,g10.3)')&
             & time_cont4_tot_max  ,  time_cont4_tot   /dble(nnod),time_cont4_tot_min  ,  time_cont4_tot   / tot_intloop
          write(*,'("INTEGRAL time in lsmpi_win_unlock ",g10.3,g10.3,g10.3,g10.3)')&
             & unlock_time_max, unlock_time/dble(nnod),unlock_time_min,          unlock_time   / tot_intloop
          write(*,'("INTEGRAL time in lsmpi_wait       ",g10.3,g10.3,g10.3,g10.3)')&
             & waiting_time_max,   waiting_time  /dble(nnod),waiting_time_min,   waiting_time  / tot_intloop
          write(*,'("INTEGRAL time in lsmpi_win_flush  ",g10.3,g10.3,g10.3,g10.3)')&
             & flushing_time_max,  flushing_time /dble(nnod),flushing_time_min,  flushing_time / tot_intloop
          write(*,'("INTEGRAL time WORK                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_w_max,phase_cntrs(PHASE_WORK_IDX)/dble(nnod),time_w_min,phase_cntrs(PHASE_WORK_IDX)/tot_intloop
          write(*,'("INTEGRAL time COMM                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_c_max,phase_cntrs(PHASE_COMM_IDX)/dble(nnod),time_c_min,phase_cntrs(PHASE_COMM_IDX)/tot_intloop
          write(*,'("INTEGRAL time IDLE                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_i_max,phase_cntrs(PHASE_IDLE_IDX)/dble(nnod),time_i_min,phase_cntrs(PHASE_IDLE_IDX)/tot_intloop
       endif
    endif



    if(DECinfo%PL>2)then
       call print_norm(integral," NORM of the integral :",print_on_rank=0)
    endif

    if(.not.local)then
       integral%access_type = AT_MASTER_ACCESS
       trafo1%access_type = AT_MASTER_ACCESS
       trafo2%access_type = AT_MASTER_ACCESS
       trafo3%access_type = AT_MASTER_ACCESS
       trafo4%access_type = AT_MASTER_ACCESS
    endif



  end subroutine get_mo_integral_par


  subroutine get_max_batch_and_scheme_ccintegral(maxsize,MinAObatch,s,MaxADA,MaxADG,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,&
        &nbuffs,nbu,set,MemToUse,use_bg_buf)
     implicit none
     type(lssetting),intent(inout) :: set
     integer, intent(in)  :: MinAObatch,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs
     integer(kind=8), intent(in) :: nbu
     real(realk), intent(in) :: MemToUse
     integer, intent(out) :: s,MaxADA,MaxADG
     integer, intent(inout) :: nbuffs
     integer(kind=long), intent(out) :: maxsize
     logical,intent(in) :: use_bg_buf
     integer :: nba, nbg, inc, i, b, e, k, magic
     integer(kind=long) :: w1size,w2size
     integer(kind=ls_mpik) :: nnod
     logical :: check_next
     nnod  = 1
     magic = 1
#ifdef VAR_MPI
     nnod  = infpar%lg_nodtot
#endif


     maxsize = 0.0E0_realk

     !get minimum memory requrements for fastest scheme
     s   = 0
     nba = MinAObatch
     nbg = MinAObatch
     inc = 1

     ! schemes 0 and 1 will return the same here!!
     w1size = get_work_array_size(1,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,nba,nbg,s,nbuffs,set)
     w2size = get_work_array_size(2,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,nba,nbg,s,nbuffs,set)

     maxsize = w1size + w2size + (i8*n1*n2)*n3*n4

     write (*,'("INFO(get_mo_integral_par): minimal memory requirements for s=0: ",g9.2," GB")')&
        &dble(maxsize*8.0E0_realk)/(1024.0**3)

     check_next = dble(maxsize*8.0E0_realk)/(1024.0**3) > MemToUse .or.&
        & (nbu < maxsize.and. use_bg_buf) .or. &
        & DECinfo%test_fully_distributed_integrals

     if(check_next)then

        s = 1
        maxsize = w1size + w2size
        write (*,'("INFO(get_mo_integral_par): minimal memory requirements for s=1: ",g9.2," GB")')&
           &dble(maxsize*8.0E0_realk)/(1024.0**3)

        check_next = dble(maxsize*8.0E0_realk)/(1024.0**3) > MemToUse .or.&
           & (nbu < maxsize.and.use_bg_buf) .or. &
           & DECinfo%test_fully_distributed_integrals

        if( check_next )then

           s = 2
           w1size = get_work_array_size(1,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,nba,nbg,s,nbuffs,set)
           w2size = get_work_array_size(2,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,nba,nbg,s,nbuffs,set)

           maxsize = w1size + w2size + (nbuffs*i8*n1s*n2s)*n3s*n4s

           write (*,'("INFO(get_mo_integral_par): minimal memory requirements for s=2: ",g9.2," GB")')&
              &dble(maxsize*8.0E0_realk)/(1024.0**3)

           check_next = dble(maxsize*8.0E0_realk)/(1024.0**3) > MemToUse .or.&
              & (nbu < maxsize.and.use_bg_buf) .or. &
              & DECinfo%test_fully_distributed_integrals

           if(check_next)then
              s      = 3
              inc    = nb/4
              nbuffs = get_nbuffs_scheme_0()
           endif

        endif
     endif


     !set requested batch sizes to the largest possible and overwrite these
     !values if this is not possible
     nba                    = nb
     nbg                    = nb

     alp: do i = MinAObatch, nb, inc
        if(s==3)then
           b = max(i/2,MinAObatch)
           e = b
        else
           b = max(i/2,MinAObatch)
           e = min(i,nb)
        endif

        do k = b, e

           w1size = get_work_array_size(1,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,i,k,s,nbuffs,set)
           w2size = get_work_array_size(2,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,i,k,s,nbuffs,set)

           maxsize = w1size + w2size

           if(s==0) maxsize = maxsize + (i8*n1*n2)*n3*n4
           if(s==2) maxsize = maxsize + (nbuffs*i8*n1s*n2s)*n3s*n4s

           if(dble(maxsize*8.0E0_realk)/(1024.0E0_realk**3) > MemToUse .or.  maxsize > nbu)then

              if(s==3)then
                 nba = max(i-inc,MinAObatch)
                 nbg = max(nba/2,MinAObatch)
              else
                 nba = max(i,MinAObatch)
                 nbg = max(min(k-1,nb),MinAObatch)
              endif

              exit alp

           endif
        enddo

     enddo alp


     if(DECinfo%manual_batchsizes)then
        nbg = min(max(DECinfo%ccsdGbatch,MinAObatch),nb)
        nba = min(max(DECinfo%ccsdAbatch,MinAObatch),nb)
     else
        !split, such that there are enough jobs for each node
        if(s/=3)then
           if((nb/nba)*(nb/nbg)<magic*nnod.and.(nba>MinAObatch).and.nnod>1)then
              do while((nb/nba)*(nb/nbg)<magic*nnod)
                 nba = max(nba - 1,MinAObatch)
                 nbg = min(max(nba/2,MinAObatch),nb)
                 if( nba == MinAObatch ) exit
              enddo
           endif

        endif

     endif


     w1size=get_work_array_size(1,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,nba,nbg,s,nbuffs,set)
     w2size=get_work_array_size(2,n1,n2,n3,n4,n1s,n2s,n3s,n4s,nb,bs,nba,nbg,s,nbuffs,set)


     maxsize = w1size + w2size

     if(s==0) maxsize = maxsize + (i8*n1*n2)*n3*n4
     if(s==2) maxsize = maxsize + (nbuffs*i8*n1s*n2s)*n3s*n4s

     print *,"REQUESTING scheme",s,"with a,g",nba,nbg,"maxsize,nbu",maxsize,nbu

     if(float(maxsize*8)/(1024.0**3) > MemToUse)then
        print*,"ERROR(get_mo_integral_par):requesting ",float(maxsize*8)/(1024.0**3),"GB available ",MemToUse,"GB"
        call lsquit("ERROR(get_mo_integral_par): the memory adaption is invalid, should not happen",-1)
     endif

     IF(use_bg_buf)THEN
        if(maxsize > nbu)then
           print*,'nbu',nbu
           print*,'maxsize',maxsize
           call lsquit("ERROR(get_mo_integral_par): the memory adaption to the bg_buffer is invalid",-1)
        endif
     ENDIF
     MaxADG = nbg
     MaxADA = nba
  end subroutine get_max_batch_and_scheme_ccintegral

  function get_work_array_size(which_array,mo1,mo2,mo3,mo4,mo1s,mo2s,mo3s,mo4s,nb,bsplit,nba,nbg,&
        &scheme,nbuf,setting) result(s)
     implicit none
     integer, intent(in) :: which_array,mo1,mo2,mo3,mo4,mo1s,mo2s,mo3s,mo4s,nb,bsplit,nba,nbg,scheme,nbuf
     type(lssetting),intent(inout) :: setting
     integer(kind=long)  :: s,maxbuf,MAX_INTEGRAL_BUF
     logical :: cd,ms
     integer :: ab,gb,starts(4),ntpm(4)
     Character :: intspec(5)

     intspec = ['R','R','R','R','C']

     select case(scheme)
     case(0)
        ms = .false.
        cd = .false.
     case(1)
        ms = .false.
        cd = .false.
     case(2)
        ms = .true.
        cd = .false.
     case(3)
        ms = .false.
        cd = .true.
     end select

     select case(which_array)
     case(1)
        if(cd)then

           maxbuf = nbuf * max(max(max(bsplit**3*mo1s,bsplit**2*mo1s*mo2s),&
              &i8*bsplit*mo1s*mo2s*mo3s),i8*mo1s*mo2s*mo3s*mo4s)

           MAX_INTEGRAL_BUF = 0

           call simulate_intloop_and_get_worksize(MAX_INTEGRAL_BUF,nb,nbg,nba,bsplit,intspec,setting)
           s = max(max(max(i8*nba*mo3,i8*nbg*mo4),maxbuf),MAX_INTEGRAL_BUF)
        else
           if(ms)then
              s = (i8*nb**2)*nba*nbg
           else
              s = max(max((i8*nb**2)*nba*nbg,(i8*mo1*mo2)*nba*nbg),(i8*mo1*mo2)*mo3*mo4)
           endif
        endif
     case(2)
        if(cd)then
           s = 1
        else
           if(ms)then
              s = (i8*mo1s*nb)*nba*nbg+(i8*mo1s*mo2s)*nba*nbg+(i8*mo1s*mo2s)*mo3s*nbg
           else
              s = max(max((i8*nb**2)*nba*nbg,(i8*mo1*mo2)*nba*nbg),(i8*mo1*mo2)*mo3*mo4)
           endif
        endif
     case default
        call lsquit("ERROR(get_work_array_size):wrong selection of array",-1)
     end select
  end function get_work_array_size

  subroutine copy_stripe_from_full_matrix(Mi,Mo,f,l,d1,d2)
     implicit none
     integer, intent(in) :: f,l,d1,d2
     real(realk), intent(in) :: Mi(d1,d2)
     real(realk), intent(out) :: Mo(l, d2)

     !$OMP WORKSHARE
     Mo(:,:) = Mi(f:f+l-1,:)
     !$OMP END WORKSHARE

  end subroutine copy_stripe_from_full_matrix

end module ccintegrals

#ifdef VAR_MPI
#ifdef MOD_UNRELEASED
!> Purpose: Intermediate routine for the slaves, they get data
!           from the local master and then call the routine to 
!           calculate MO integrals (non-T1 transformed)
!
!> Author:  Pablo Baudin
!> Date:    December 2013
subroutine cc_gmo_data_slave()

  use dec_typedef_module
  use ccintegrals
  use daltoninfo
  use typedeftype, only: lsitem
  use decmpi_module, only: mpi_communicate_get_gmo_data

  implicit none

  !> number of orbitals:
  integer :: nb, no, nv
  !> number of MO batch
  integer :: nbatch
  !> SCF transformation matrices:
  real(realk), pointer  :: Co(:,:), Cv(:,:)
  !> performed MO-based CCSD calculation ?
  logical :: mo_ccsd
  !> array with gmo on output:
  type(tensor) :: pgmo_diag, pgmo_up
  !> variables used for MO batch and integral transformation
  type(MObatchInfo) :: MOinfo
  !> LS item information
  type(lsitem) :: MyLsItem


  call mpi_communicate_get_gmo_data(mo_ccsd,MyLsItem,Co,Cv, &
       & pgmo_diag,pgmo_up,nb,no,nv,nbatch)

  ! the slave call the routine to get MO int.
  call get_t1_free_gmo(mo_ccsd,MyLsItem,Co,Cv, &
       & pgmo_diag,pgmo_up,nb,no,nv,MOinfo)

  ! deallocate slave stuff:
  call mem_dealloc(Co)
  call mem_dealloc(Cv)
  call ls_free(MyLsItem)
  call mem_dealloc(MOinfo%dimInd1)
  call mem_dealloc(MOinfo%dimInd2)
  call mem_dealloc(MOinfo%StartInd1)
  call mem_dealloc(MOinfo%StartInd2)
  call mem_dealloc(MOinfo%dimTot)
  call mem_dealloc(MOinfo%tileInd)

end subroutine cc_gmo_data_slave
#endif

subroutine get_mo_integral_par_slave()
  use dec_typedef_module
  use daltoninfo
  use tensor_type_def_module
  use typedeftype, only: lsitem
  use decmpi_module, only: wake_slaves_for_simple_mo
  use ccintegrals, only : get_mo_integral_par

  implicit none
  type(tensor) :: integral,trafo1,trafo2,trafo3,trafo4
  type(lsitem) :: mylsitem
  logical :: c
  character :: is(5)


  call wake_slaves_for_simple_mo(integral,trafo1,trafo2,trafo3,trafo4,mylsitem,c)

  ! Set integral info
  ! *****************
  !R = Regular Basis set on the 1th center 
  !R = Regular Basis set on the 2th center 
  !R = Regular Basis set on the 3th center 
  !R = Regular Basis set on the 4th center 
  !C = Coulomb operator
  !E = Long-Range Erf operator
  if (mylsitem%setting%scheme%CAM) then
     is = ['R','R','R','R','E'] 
  else
     is = ['R','R','R','R','C'] 
  endif

  call get_mo_integral_par(integral,trafo1,trafo2,trafo3,trafo4,mylsitem,is,.false.,c)
  call ls_free(mylsitem)

end subroutine get_mo_integral_par_slave
#endif
