!> @file
!> Simple integrals related
!> \author Marcin Ziolkowski, (Pablo Baudin: added subroutines to get 
!          (non-T1-transformed MO int. used for MO-based CCSD and RPA)
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
#ifdef VAR_ICHOR
   use IchorErimoduleHost
#endif

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
#ifdef VAR_MPI
  use decmpi_module
#endif

  interface getL
     module procedure getL_simple
     module procedure getL_simple_from_gmo
     module procedure getL_diff
  end interface

  private :: get_mem_t1_free_gmo, get_mem_MO_CCSD_residual, &
       & get_AO_batches_size_rpa, get_mem_gmo_RPA, gao_to_gmo, gao_to_govov, &
       & get_MO_batches_info, pack_and_add_gmo

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

    write(DECinfo%output,'(a)') 'info :: calculating two-electron integrals'
    g_ao = array4_init_standard(ao_dims)
    ! KK Quick fix: Filename associated with g_ao
    g_ao%filename = 'gao'
#ifdef VAR_ICHOR
    SameMOL = .TRUE.
    iprint = 0
    call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,&
         & INTSPEC,SameMOL)
    call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,&
         & nbasis,nbasis,nbasis,nbasis,g_ao%val,INTSPEC,.TRUE.,&
         & 1,1,1,1,1,1,1,1,.FALSE.,nbasis,nbasis,nbasis,nbasis,.FALSE.)
    call FREE_SCREEN_ICHORERI()
#else
    call ii_get_4center_eri(DECinfo%output,DECinfo%output,&
         & mylsitem%setting,g_ao%val,nbasis,nbasis,nbasis,nbasis,intspec)
#endif
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
#ifdef VAR_ICHOR
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
         & MoTrans,ndim(1),ndim(2),ndim(3),ndim(4),NoSymmetry)
    !Free screening info
    call FREE_SCREEN_ICHORERI()
#else
    !Use Thermite code to calculate Integrals     
    ! Set integral screening
    doscreen = mylsitem%setting%scheme%cs_screen .OR. mylsitem%setting%scheme%ps_screen
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mylsitem%setting,1,1,intspecConvert)
    IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
    IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabRHS

    ! Get AO integrals
    call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output,Mylsitem%SETTING,&
         & gao,1,1,ndim(3),ndim(4),ndim(1),ndim(2),ndim(3),ndim(4),.true.,intSpecConvert)
    call free_decscreen(DECSCREEN)
    nullify(mylsitem%setting%LST_GAB_RHS)
    nullify(mylsitem%setting%LST_GAB_LHS)
#endif

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
    logical :: doscreen
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
    CALL II_get_h1_mixed(DECinfo%output,DECinfo%output,MyLsitem%SETTING,h,AO(1),AO(2))
    CALL II_get_coulomb_mat_mixed(DECinfo%output,DECinfo%output,MyLsitem%SETTING,(/D/),Jarr,ndmat,&
         &                            AO(1),AO(2),AO(3),AO(4),Oper)

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
    logical :: Dsym
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
    CALL ii_get_exchange_mat_mixed(DECinfo%output,DECinfo%output,MyLsitem%SETTING,(/D/),ndmat,Dsym,&
         &                             Karr,AO(1),AO(3),AO(2),AO(4),Oper)
    call mat_assign(K,Karr(1)) 
    call mat_free(Karr(1))

  end subroutine get_AO_K


#ifdef MOD_UNRELEASED
  !> Purpose: calculate AO int. in batches and transform them to
  !           full MO basis (non T1-transformed)
  !           The batches are then packed using permutational
  !           symmetry and are kept in memory (PDM if MPI)
  !           If the routine is call for RPA then only govov is 
  !           calculated without batching and packing.
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine get_t1_free_gmo(mo_ccsd,mylsitem,Co,Cv,govov,pgmo_diag,pgmo_up, &
       & nb,no,nv,CCmodel,MOinfo)

    implicit none

    !> number of orbitals:
    integer, intent(in) :: nb, no, nv, CCmodel
    !> SCF transformation matrices:
    real(realk), pointer, intent(in) :: Co(:,:), Cv(:,:)
    !> performed MO-based CCSD calculation ?
    logical, intent(inout) :: mo_ccsd
    !> array with packed gmo on output:
    ! (intent in needed for the slaves)
    type(tensor), intent(inout) :: govov
    type(tensor), intent(inout) :: pgmo_diag, pgmo_up

    !> variables used for MO batch and integral transformation
    integer :: ntot ! total number of MO
    real(realk), pointer :: Cov(:,:), CP(:,:), CQ(:,:)
    real(realk), pointer :: gmo(:), tmp1(:), tmp2(:)
    integer(kind=long)   :: gmosize, tmp_size
    integer :: Nbatch, PQ_batch, dimP, dimQ, idb, iub
    integer :: P_sta, P_end, Q_sta, Q_end
    type(MObatchInfo), intent(out) :: MOinfo
    logical :: local_moccsd

    !> variables used for AO batch construction and AO integral calculation
    real(realk), pointer :: gao(:)
    integer(kind=long)   :: gaosize
    integer :: alphaB, gammaB, dimAlpha, dimGamma
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
    integer :: iorb, idx, K
#ifdef VAR_ICHOR
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
#else
    integer, pointer :: batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchAlpha(:),orb2batchGamma(:)
    integer, pointer :: batchsizeGamma(:), batchindexGamma(:)
    ! Screening integrals stuff:
    type(DECscreenITEM) :: DecScreen
#endif
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
    logical :: master, local, gdi_lk, gup_lk
    integer :: myload, win
    integer(kind=ls_mpik) :: ierr, myrank, nnod, dest
    integer, pointer      :: tasks(:)

    !> LS item with information needed for integrals
    type(lsitem), intent(inout) :: MyLsItem


    call time_start_phase(PHASE_WORK)
    ntot = no + nv

    ! Set integral info
    ! *****************
    INTSPEC(1)  = 'R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)  = 'R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)  = 'R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)  = 'R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)  = 'C' !C = Coulomb operator
#ifdef VAR_ICHOR
    iprint = 0           !print level for Ichor Integral code
    MoTrans = .FALSE.    !Do not transform to MO basis! 
    NoSymmetry = .FALSE. !Use Permutational Symmetry! 
    SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
    !Determine the full number of AO batches - not to be confused with the batches of AOs
    !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
    iAO = 1
    call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
#else
    doscreen = MyLsItem%setting%scheme%cs_screen.OR. &
         & MyLsItem%setting%scheme%ps_screen
#endif
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
#ifdef VAR_ICHOR
    nullify(AOGammabatchinfo)
    nullify(AOalphabatchinfo)    
#else
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
#endif

    nullify(Cov)
    nullify(CP)
    nullify(CQ)
    nullify(gao)
    nullify(gmo)
    nullify(tmp1)
    nullify(tmp2)
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
       select case(CCmodel)

       case(MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT)
          call get_MO_and_AO_batches_size(mo_ccsd,local_moccsd,ntot,nb,no,nv, &
               & dimP,Nbatch,MaxAllowedDimAlpha,MaxAllowedDimGamma,MyLsItem,.false.)
          if (.not.mo_ccsd) return

          if (print_debug) then
             if (local_moccsd) then 
                write(DECinfo%output,*) 'MO-CCSD: local scheme'
             else if (.not.local) then
                write(DECinfo%output,*) 'MO-CCSD: PDM scheme'
             else
                write(DECinfo%output,*) 'MO-CCSD: non-MPI scheme'
             end if
             write(DECinfo%output,'(a,I4,a,I4)') ' BATCH: Number of MO batches      = ', &
                  & Nbatch*(Nbatch+1)/2, ' with maximum size', dimP
          end if

          ! Initialize gmo arrays:
          call init_gmo_arrays(ntot,dimP,Nbatch,local,local_moccsd,pgmo_diag,pgmo_up)

       case(MODEL_RPA)
          call get_AO_batches_size_rpa(ntot,nb,no,nv,MaxAllowedDimAlpha, &
               & MaxAllowedDimGamma,MyLsItem)
       case default
          call lsquit('only RPA, CCSD and CCSD(T) model should use this routine',DECinfo%output)
       end select

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
            & pgmo_diag,pgmo_up,nb,no,nv,Nbatch,ccmodel)
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
#ifdef VAR_ICHOR
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
#else
    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchGamma,nb)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma, &
         & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma, &
         & nbatchesGamma,orb2BatchGamma,'R')
#endif

    if (print_debug) write(DECinfo%output,'(a,I4,a,I4)') & 
         & ' BATCH: Number of Gamma batches   = ', nbatchesGamma, &
         & ' with maximum size', MaxActualDimGamma 


#ifndef VAR_ICHOR
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
#endif


    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

#ifdef VAR_ICHOR
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
#else
    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nb)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha, &
         & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha, &
         & nbatchesAlpha,orb2BatchAlpha,'R')
#endif

    if (print_debug) write(DECinfo%output,'(a,I4,a,I4)') & 
         & ' BATCH: Number of Alpha batches   = ', nbatchesAlpha, &
         & ' with maximum size',MaxActualDimAlpha


#ifndef VAR_ICHOR
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
#endif

    ! **************************************
    ! * Allocate Memory to working arrays  *
    ! **************************************

    ! AO integral allocation:
    gaosize = int(i8*nb*nb*MaxActualDimAlpha*MaxActualDimGamma,kind=long)
    call mem_alloc(gao,gaosize)


    if(ccmodel == MODEL_RPA) then
       ! working arrays
       gmosize = int(i8*no*nv*no*nv,kind=long)
       call mem_alloc(gmo,gmosize)
       gmo = 0.0_realk
       tmp_size = max(nb*MaxActualDimAlpha*MaxActualDimGamma, MaxActualDimGamma*no*nv)
       tmp_size = int(i8*tmp_size*no, kind=long)
       call mem_alloc(tmp1, tmp_size)
       tmp_size = int(i8*MaxActualDimAlpha*MaxActualDimGamma*no*nv, kind=long)
       call mem_alloc(tmp2, tmp_size)
    else
       ! Get full MO coeficients:
       call mem_alloc(Cov,nb,ntot)
       Cov(:,:no)       = Co
       Cov(:,no+1:ntot) = Cv

       ! gmo batch array:
       gmosize = int(i8*dimP*dimP*ntot*ntot,kind=long)
       call get_MO_batches_info(MOinfo, dimP, ntot, Nbatch)
       call mem_alloc(CP,MaxActualDimAlpha,dimP)
       call mem_alloc(CQ,MaxActualDimGamma,dimP)
       call mem_alloc(gmo,gmosize)

       ! working arrays
       tmp_size = max(nb*MaxActualDimAlpha*MaxActualDimGamma, ntot*MaxActualDimGamma*dimP)
       tmp_size = int(i8*ntot*tmp_size, kind=long)
       call mem_alloc(tmp1, tmp_size)
       tmp_size = max(MaxActualDimAlpha*MaxActualDimGamma, dimP*dimP)
       tmp_size = int(i8*ntot*ntot*tmp_size, kind=long)
       call mem_alloc(tmp2, tmp_size)
    endif


    ! Sanity checks for matrix sizes which need to be filled:
    if (gaosize>MaxInt) then
       call lsquit("ERROR(CCSD):matrix sizes too large, &
            & please recompile with 64bit integers",-1)
    endif


#ifdef VAR_ICHOR
    !Calculate Screening integrals 
    SameMOL = .TRUE. !Specifies same molecule on all centers 
    call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
#else
    ! *******************************************************
    ! *  This subroutine builds the full screening matrix.
    call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting, &
         & nbatchesAlpha,nbatchesGamma,INTSPEC)

    if (mylsitem%setting%scheme%cs_screen .or. mylsitem%setting%scheme%ps_screen) then
       call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
            & nb,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
    end if
    ! *******************************************************
#endif

#ifdef VAR_MPI
    ! Calculate the batches for a good load balance
    call mem_alloc(tasks,nbatchesAlpha*nbatchesGamma)

    myload = 0
    tasks  = 0
#ifdef VAR_ICHOR
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
#else
    call distribute_mpi_jobs(tasks,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
         &batchdimGamma,myload,nnod,myrank,4,no,nv,nb,batch2orbAlpha,&
         &batch2orbGamma)
#endif
#endif
    myload = 0

    fullRHS = (nbatchesGamma.eq.1).and.(nbatchesAlpha.eq.1)

    !**********************************
    ! Begin the loop over gamma batches
    !**********************************


    BatchGamma: do gammaB = 1,nbatchesGamma            ! batches of AO batches
#ifdef VAR_ICHOR
       dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
       GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
       GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
       AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch needed by MAIN_ICHORERI_DRIVER
       AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch needed by MAIN_ICHORERI_DRIVER
#else

       dimGamma   = batchdimGamma(gammaB)                         ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd   = batch2orbGamma(gammaB)%orbindex(dimGamma)     ! Last index in gamma batch
#endif
       !**********************************
       ! Begin the loop over alpha batches
       !**********************************

       BatchAlpha: do alphaB = 1,nbatchesAlpha         ! batches of AO batches

          if (nnod>1) then 
             ! check if the current job is to be done by current node
             if (tasks(alphaB + (gammaB-1)*nbatchesAlpha)/=myrank) cycle
          end if

#ifdef VAR_ICHOR
          dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
          AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
          AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
          AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
          AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
#else
          dimAlpha   = batchdimAlpha(alphaB)                        ! Dimension of alpha batch
          AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)           ! First index in alpha batch
          AlphaEnd   = batch2orbAlpha(alphaB)%orbindex(dimAlpha)    ! Last index in alpha batch
#endif
          myload     = myload + dimAlpha*dimGamma


#ifdef VAR_ICHOR
          call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nb,nb,dimAlpha,dimGamma,&
               & gao,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
               & AOGammaStart,AOGammaEnd,MoTrans,nb,nb,dimAlpha,dimGamma,NoSymmetry)
#else
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
               & dimGamma,fullRHS,INTSPEC)
#endif
          if (ccmodel == MODEL_RPA) then


             !call gao_to_govov(govov%elm1,gao,Co,Cv,nb,no,nv,AlphaStart,dimAlpha, &
             !     & GammaStart,dimGamma,tmp1,tmp2)
             call gao_to_govov(gmo,gao,Co,Cv,nb,no,nv,AlphaStart,dimAlpha, &
                  & GammaStart,dimGamma,tmp1,tmp2)

          else
             idb = 0
             iub = 0
             ! Loop over MO batches:
             BatchPQ: do PQ_batch = 1, MOinfo%nbatch


                P_sta  = MOinfo%StartInd1(PQ_batch)
                dimP   = MOinfo%DimInd1(PQ_batch)
                Q_sta  = MOinfo%StartInd2(PQ_batch)
                dimQ   = MOinfo%DimInd2(PQ_batch)

                call gao_to_gmo(gmo,gao,Cov,CP,CQ,nb,ntot,AlphaStart,dimAlpha, &
                     & GammaStart,dimGamma,P_sta,dimP,Q_sta,dimQ,tmp1, &
                     & tmp2,pgmo_diag,pgmo_up,gdi_lk,gup_lk,win,dest)

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
                   call pack_and_add_gmo(gmo,pgmo_diag,idb,ntot,dimP,dimQ,.true.,tmp2)
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
                   call pack_and_add_gmo(gmo,pgmo_up,iub,ntot,dimP,dimQ,.false.,tmp2)
                end if

             end do BatchPQ

          end if


       end do BatchAlpha
    end do BatchGamma

    ! Free integral stuff
    ! *******************
#ifdef VAR_ICHOR
    call FREE_SCREEN_ICHORERI()
    call mem_dealloc(AOGammabatchinfo)
    call mem_dealloc(AOAlphabatchinfo)
#else
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
#endif

    if (ccmodel==MODEL_RPA) then 
       !call tensor_scatter(1.0E0_realk,gmo,0.0E0_realk,govov,i8*no*nv*no*nv)
       if(master) then
          !  call print_norm(gmo,i8*no*no*nv*nv)
          call tensor_convert(gmo,govov)
          !  call print_norm(govov)
       endif
       !call daxpy(ncopy,1.0E0_realk,gmo,1,govov%elm1,1)
    endif

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
    if (.not.local_moccsd.and.ccmodel/=MODEL_RPA) then
       call time_start_phase(PHASE_IDLE)
       call lsmpi_barrier(infpar%lg_comm)
       call time_start_phase(PHASE_WORK)
    end if
#endif

    ! Free matrices:
    call mem_dealloc(gao)
    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
    call mem_dealloc(gmo)
    if (ccmodel/=MODEL_RPA) then 
       call mem_dealloc(Cov)
       call mem_dealloc(CP)
       call mem_dealloc(CQ)
    end if

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

    real(realk) :: MemNeed, MemFree
    integer(kind=long) :: min_mem
    integer :: MinAOBatch, MinMOBatch, na, ng, nnod, magic,iAO

    MinMOBatch = min(15,ntot)
    dimMO = MinMOBatch
    local = .false.
    nnod  = 1
#ifdef VAR_MPI
    nnod  = infpar%lg_nodtot
#endif

    !===========================================================
    ! Get MO batch size depending on MO-ccsd residual routine.
    call get_currently_available_memory(MemFree)
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
          call get_mem_MO_CCSD_residual(local,MemNeed,ntot,nb,no,nv,dimMO)

          ! if not enough mem then switch to full PDM scheme:
          if (MeMNeed>0.8E0_realk*MemFree) then
             local = .false.
          end if
       end if
    end if

    call get_mem_MO_CCSD_residual(local,MemNeed,ntot,nb,no,nv,dimMO)

    do while ((MemNeed<0.8E0_realk*MemFree).and.(dimMO<=ntot))
       dimMO = dimMO + 1
       call get_mem_MO_CCSD_residual(local,MemNeed,ntot,nb,no,nv,dimMO)
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
    call get_mem_MO_CCSD_residual(local,MemNeed,ntot,nb,no,nv,dimMO) 
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
#ifdef VAR_ICHOR
    !Determine the minimum allowed AObatch size MinAObatch
    !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
    !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
    !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
    !'R'  !Specifies that it is the Regular AO basis that should be batched
    iAO = 1 !the center that the batching should occur on.  
    call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
#else
    call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
#endif
    call get_mem_t1_free_gmo(MemNeed,ntot,nb,no,nv,dimMO,Nbatch, &
         & MinAObatch,MinAObatch,MinAObatch)

    MaxGamma = MinAObatch
    MaxAlpha = MinAObatch
    do while ((MemNeed<0.8E0_realk*MemFree).and.(MaxGamma<=nb)) 
       MaxGamma = MaxGamma + 1
       call get_mem_t1_free_gmo(MemNeed,ntot,nb,no,nv,dimMO,Nbatch, &
            & MaxAlpha,MaxGamma,MinAObatch)  
    end do
    if (MaxGamma>=nb) then
       MaxGamma = nb
    else if (MaxGamma<=MinAObatch) then
       MaxGamma = MinAObatch
    else 
       MaxGamma = MaxGamma - 1
    end if
    do while ((MemNeed<0.8E0_realk*MemFree).and.(MaxAlpha<=nb)) 
       MaxAlpha = MaxAlpha + 1
       call get_mem_t1_free_gmo(MemNeed,ntot,nb,no,nv,dimMO,Nbatch, &
            & MaxAlpha,MaxGamma,MinAObatch)  
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
    call get_mem_t1_free_gmo(MemNeed,ntot,nb,no,nv,dimMO,Nbatch, &
         & MaxAlpha,MaxGamma,MinAObatch)  
    if ((MemFree-MemNeed)<=0.0E0_realk) then
       mo_ccsd = .false.
       write(DECinfo%output,*) 'WARNING: Insufficient memory in MO-based CCSD, &
            & back to standard algorithm.'
       write(DECinfo%output,'(a,F12.5,a)') '   Available memory:',MemFree,' GB'
       write(DECinfo%output,'(a,F12.5,a)') '   Required memory :',MemNeed,' GB'
       return
    end if

  end subroutine get_MO_and_AO_batches_size


  !> Purpose: Get memory required in get_t1_free_gmo depending on AO 
  !           batch dimension.
  !           Get min. required memory: AlphaDim = GammaDim = MinDimAO
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  subroutine get_mem_t1_free_gmo(MemOut,M,N,O,V,X,nMOB,AlphaDim,GammaDim,MinDimAO)

    implicit none 

    ! M: tot number of MO
    ! N: tot number of AO
    ! O: number of occ. orbs.
    ! V: number of virt. orbs.
    ! X: dimension of MO batch.
    ! nMOB: number of MO batches.
    integer,  intent(in) :: M, N, O, V, X, nMOB
    !> AO stuff:
    integer, intent(in) :: AlphaDim, GammaDim, MinDimAO
    !> memory needed:
    real(realk), intent(inout) :: MemOut
    ! intermediate memory:
    integer :: nnod
    integer(kind=long) :: MemNeed, nTileMax

    nnod = 1
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
#endif

    ! Transfo. matrices:
    MemNeed = N*M + AlphaDim*X + GammaDim*X

    ! AO stuff:
    MemNeed = MemNeed + 4*N + N*N*AlphaDim*GammaDim 

    ! Packed gmo diag blocks:
    nTileMax = (nMOB-1)/nnod + 3
    MemNeed = MEmNeed + nTileMax*X*(X+1)*M*(M+1)/4
    ! Packed gmo upper blocks:
    nTileMax = (nMOB*(nMOB-1)/2 - 1)/nnod + 3
    MemNeed = MEmNeed + nTileMax*X*X*M*(M+1)/2

    ! MO stuff:
    MemNeed = MemNeed + X*X*M*M + 5*nMOB*nMOB + 1

    ! Working arrays:
    MemNeed = MemNeed + max(M*N*AlphaDim*GammaDim, M*M*GammaDim*X)
    MemNeed = MemNeed + max(M*M*AlphaDim*GammaDim, M*M*X*X)

    MemOut = MemNeed*8.0E0_realk/(1.024E3_realk**3) 

  end subroutine get_mem_t1_free_gmo


  !> Purpose: Get memory required in get_ccsd_residual_mo_ccsd 
  !           depending on MO batch dimension.
  !           Get min. required memory when X = 1 
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  subroutine get_mem_MO_CCSD_residual(local,MemOut,M,N,O,V,X)

    implicit none 

    ! M: tot number of MO
    ! N: tot number of AO
    ! O: number of occ. orbs.
    ! V: number of virt. orbs.
    ! X: dimension of MO batch.
    integer,  intent(in) :: M, N, O, V, X
    !> use local scheme?
    logical :: local
    !> memory needed:
    real(realk), intent(out) :: MemOut
    !> intermediate memory:
    integer :: nnod, nMOB
    integer(kind=long) :: nTileMax, MemNeed

    nMOB = (M-1)/X + 1
    nnod = 1
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
#endif
    if (local) nnod = 1

    ! Packed gmo diag blocks:
    nTileMax = (nMOB-1)/nnod + 3
    MemNeed = nTileMax*X*(X+1)*M*(M+1)/4

    ! Packed gmo upper blocks:
    nTileMax = (nMOB*(nMOB-1)/2 - 1)/nnod + 3
    MemNeed = MemNeed + nTileMax*X*X*M*(M+1)/2

    ! Working arrays:
    MemNeed = MemNeed + max(O**4, V*O**3, V*V*O*O, X*X*M*M, X*O*O*V, X*O*V*V)
    MemNeed = MemNeed + max(X*X*M*M, O*O*V*M, O*O*X*M, O**4)
    MemNeed = MemNeed + max(X*O*V*M, O*O*V*V, X*X*M*M, X*O*O*M)

    ! T1-Transfo. matrices:
    MemNeed = MemNeed + V*M + O*M

    ! Batch of MO int:
    MemNeed = MemNeed + X*X*M*M

    ! Intermediates (B2prep, u2, G_Pi, H_aQ):
    MemNeed = MemNeed + O**4 + O*O*V*V + X*O + V*X

    ! T1-transformed integrals:
    MemNeed = MemNeed + O**4 + 2*V*O**3 + 3*O*O*V*V 

    ! Fock Matrix:
    MemNeed = MemNeed + 3*N*N

    MemOut = MemNeed*8.0E0_realk/(1.024E3_realk**3) 

  end subroutine get_mem_MO_CCSD_residual


  !> Purpose: Calculate Max. AO batches based on Min. AO batches size, 
  !           available memory and mem. requirements in get_t1_free_gmo routine.
  !
  !> Author:  Johannes Rekkedal
  !> Date:    January 2014
  subroutine get_AO_batches_size_rpa(ntot,nbas,nocc,nvir,MaxAlpha, &
       & MaxGamma,MyLsItem)

    implicit none

    !> number of orbitals:
    integer, intent(in) :: ntot, nbas, nocc, nvir
    !> AO batches stuff:
    integer, intent (inout) :: MaxAlpha, MaxGamma
    type(lsitem), intent(inout) :: MyLsItem

    real(realk) :: MemNeed, MemFree
    integer(kind=long) :: min_mem
    integer :: MinAOBatch

    ! Get minimum mem. required to get gmo:
    call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
    call get_currently_available_memory(MemFree)
    call get_mem_gmo_RPA(MemNeed,ntot,nbas,nocc,nvir,MinAObatch,MinAObatch,MinAObatch)


    MaxGamma = MinAObatch
    MaxAlpha = MinAObatch
    do while ((MemNeed<0.8E0_realk*MemFree).and.(MaxGamma<=nbas)) 
       MaxGamma = MaxGamma + 1
       call get_mem_gmo_RPA(MemNeed,ntot,nbas,nocc,nvir,MaxAlpha,MaxGamma,MinAObatch)  
    end do

    if (MaxGamma>=nbas) then
       MaxGamma = nbas
    else if (MaxGamma<=MinAObatch) then
       MaxGamma = MinAObatch
    else 
       MaxGamma = MaxGamma - 1
    end if

    do while ((MemNeed<0.8E0_realk*MemFree).and.(MaxAlpha<=nbas)) 
       MaxAlpha = MaxAlpha + 1
       call get_mem_gmo_RPA(MemNeed,ntot,nbas,nocc,nvir,MaxAlpha,MaxGamma,MinAObatch)  
    end do

    if (MaxAlpha>=nbas) then
       MaxAlpha = nbas
    else if (MaxAlpha<=MinAObatch) then
       MaxAlpha = MinAObatch
    else 
       MaxAlpha = MaxAlpha - 1
    end if

    ! sanity check:
    call get_mem_gmo_RPA(MemNeed,ntot,nbas,nocc,nvir,MaxAlpha,MaxGamma,MinAObatch)  
    if ((MemFree-MemNeed)<=0.0E0_realk) then
       call lsquit('Not enough memory in RPA (MO int calc.)', DECinfo%output)
    end if

  end subroutine get_AO_batches_size_RPA


  !> Purpose: Get memory required in get_t1_free_gmo for RPA model
  !
  !> Author:  Johannes Rekkedal
  !> Date:    January 2014
  subroutine get_mem_gmo_RPA(MemOut,M,N,O,V,AlphaDim,GammaDim,MinDimAO)

    implicit none 

    ! M: tot number of MO
    ! N: tot number of AO
    ! O: number of occ. orbs.
    ! V: number of virt. orbs.
    integer,  intent(in) :: M, N, O, V
    !> AO stuff:
    integer, intent(in) :: AlphaDim, GammaDim, MinDimAO
    !> memory needed:
    real(realk), intent(inout) :: MemOut
    ! intermediate memory:
    integer :: MemNeed


    ! AO stuff:
    MemNeed = 4*N + N*N*AlphaDim*GammaDim 

    ! Working arrays:
    MemNeed = MemNeed + max(N*AlphaDim*GammaDim*O, GammaDim*O*V*O)
    MemNeed = MemNeed + AlphaDim*GammaDim*O*V

    MemOut = MemNeed*8.0E0_realk/(1.024E3_realk**3) 

  end subroutine get_mem_gmo_RPA


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

    character(4) :: at
    integer :: pgmo_dims

    ! define type of array:
    if (local_moccsd) then
       at = 'RTAR'
    else 
       at = 'TDAR'
    end if

    ! Declare one array for the diagonal batches:
    pgmo_dims = ntot*(ntot+1)*dimMO*(dimMO+1)/4
    call tensor_minit(pgmo_diag,[pgmo_dims,Nbat],2,local=mpi,atype=at,tdims=[pgmo_dims,1])
    call tensor_zero(pgmo_diag)

    ! Declare one array for the upper diagonal batches if necesarry
    if (Nbat>1) then
      pgmo_dims = ntot*(ntot+1)*dimMO*dimMO/2
      call tensor_minit(pgmo_up,[pgmo_dims,Nbat*(Nbat-1)/2],2,local=mpi,atype=at, &
                & tdims=[pgmo_dims,1])
      call tensor_zero(pgmo_up)
    end if

  end subroutine init_gmo_arrays


  !> Purpose: Transform AO int. into MO int. in batches
  !           
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine gao_to_gmo(gmo,gao,Cov,CP,CQ,nb,ntot,AlphaStart,dimAlpha, &
       & GammaStart,dimGamma,P_sta,dimP,Q_sta,dimQ,tmp1,tmp2, &
       & pgmo_diag,pgmo_up,gdi_lk,gup_lk,win,dest)

    implicit none

    integer, intent(in) :: nb, AlphaStart, dimAlpha, GammaStart, dimGamma
    integer, intent(in) :: ntot, P_sta, dimP, Q_sta, dimQ
    real(realk), intent(inout) :: gmo(dimP*dimQ*ntot*ntot)
    real(realk), intent(in) :: gao(nb*nb*dimAlpha*dimGamma), Cov(nb,ntot)
    !> MPI related:
    type(tensor), intent(in) :: pgmo_diag, pgmo_up
    logical, intent(inout)  :: gdi_lk, gup_lk
    integer, intent(in) :: win
    integer(kind=ls_mpik), intent(in) :: dest

    real(realk) :: CP(dimAlpha,dimP), CQ(dimGamma,dimQ)
    real(realk) :: tmp1(:), tmp2(:)
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
         & gao,nb,Cov,nb,0.0E0_realk,tmp1,nb*dimAlpha*dimGamma)

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
         & tmp1,nb,Cov,nb,0.0E0_realk,tmp2,dimAlpha*dimGamma*ntot)

    ! transfo alphaB to P_batch => [gammaB r s, P]
    call dgemm('t','n',dimGamma*ntot*ntot,dimP,dimAlpha,1.0E0_realk, &
         & tmp2,dimAlpha,CP,dimAlpha,0.0E0_realk,tmp1,dimGamma*ntot*ntot)

    ! transfo gammaB to Q_batch => [r s P, Q]
    call dgemm('t','n',ntot*ntot*dimP,dimQ,dimGamma,1.0E0_realk, &
         & tmp1,dimGamma,CQ,dimGamma,0.0E0_realk,tmp2,ntot*ntot*dimP)

    ! transpose matrix => [P_batch, Q_batch, r, s]
    call mat_transpose(ntot*ntot,dimP*dimQ,1.0E0_realk,tmp2,0.0E0_realk,gmo)

  end subroutine gao_to_gmo

  !> Purpose: Transform AO int. into MO (occ,virt,occ,virt) in batches
  !           
  !> Author:  Johannes Rekkedal
  !> Date:    December 2013
  subroutine gao_to_govov(govov,gao,Co,Cv,nb,no,nv,A_sta,dimAlpha, &
       & G_Sta,dimGamma,tmp1,tmp2)

    implicit none

    ! array dimensions:
    integer, intent(in) :: nb, no, nv
    integer, intent(in) :: A_Sta, dimAlpha, G_Sta, dimGamma

    !> MO integral:
    real(realk), intent(inout) :: govov(no*nv*no*nv) 
    !> Batch of AO integral:
    real(realk), intent(in) :: gao(dimAlpha*nb*dimGamma*nb)
    !> Transfo. matrices:
    real(realk), intent(in) :: Co(nb,no),Cv(nb,nv)
    !> working arrays:
    real(realk) :: tmp1(:), tmp2(:)

    integer :: A_end, G_end

    A_end = A_sta +dimAlpha - 1
    G_end = G_sta +dimGamma - 1



    ! we have (beta delta alpha gamma)

    ! transfo Beta to j => [delta alphaB gammaB, j]
    call dgemm('t','n',nb*dimAlpha*dimGamma,no,nb,1.0E0_realk, &
         & gao,nb,Co,nb,0.0E0_realk,tmp1,nb*dimAlpha*dimGamma)

    ! transfo delta to b => [alphaB gammaB j, b]
    call dgemm('t','n',dimAlpha*dimGamma*no,nv,nb,1.0E0_realk, &
         & tmp1,nb,Cv,nb,0.0E0_realk,tmp2,dimAlpha*dimGamma*no)

    ! transfo alphaB to i => [gammaB j b, i]
    call dgemm('t','n',dimGamma*no*nv,no,dimAlpha,1.0E0_realk,tmp2,dimAlpha, &
         & Co(A_sta:A_end,:),dimAlpha,0.0E0_realk,tmp1,dimGamma*no*nv)

    ! transfo gammaB to a => [j b i, a]
    call dgemm('t','n',no*nv*no,nv,dimGamma,1.0E0_realk,tmp1,dimGamma, &
         & Cv(G_sta:G_end,:),dimGamma,1.0E0_realk,govov,no*nv*no)


  end subroutine gao_to_govov


  !> Purpose: Get information regarding MO batches for two 
  !           equivalent nested loops:
  !
  !> Author:  Pablo Baudin
  !> Date:    Novemeber 2013
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
  subroutine pack_and_add_gmo(gmo,pack_gmo,tile,ntot,dimP,dimQ,diag,tmp)

    implicit none

    !> array with one batch of partial MO int.:
    real(realk), intent(in) :: gmo(:)
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
    real(realk), intent(inout) :: tmp(:)

    integer :: s, r, rs, q, ibatch, ipack, nnod
    integer(kind=long) :: ncopy

    nnod = 1
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
#endif

    ! 1st case: current batch corresponds to diagonal block, we 
    !           keep only the upper triangular part of the batch.
    if (diag) then

       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ntot,gmo,tmp,dimP,dimQ)&
       !$OMP PRIVATE(q,r,s,rs,ibatch,ipack)
       do s=1,ntot
          do r=1,s
             rs = r + (s-1)*ntot
             do q=1,dimQ
                ibatch = 1 + (q-1)*dimP + (rs-1)*dimP*dimQ
                ipack = 1 + q*(q-1)/2 + (s*(s-1)/2 + r-1)*(dimQ*(dimQ+1)/2)
                call dcopy(q,gmo(ibatch),1,tmp(ipack),1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       ncopy = dimQ*(dimQ+1)/2 + (ntot*(ntot-1)/2 + ntot-1)*(dimQ*(dimQ+1)/2)
       ! accumulate tile
       if (nnod>1.and.pack_gmo%itype==TT_TILED_DIST) then
          call time_start_phase(PHASE_COMM)
          call tensor_accumulate_tile(pack_gmo,tile,tmp(1:ncopy),ncopy,lock_set=.true.)
          call time_start_phase(PHASE_WORK)
       else if (nnod>1.and.pack_gmo%itype==TT_TILED) then
          call daxpy(ncopy,1.0E0_realk,tmp,1,pack_gmo%ti(tile)%t(:),1)
       else
          call daxpy(ncopy,1.0E0_realk,tmp,1,pack_gmo%elm2(:,tile),1)
       end if

       ! 2nd case: current batch corresponds to an upper diagonal block,
       !           we keep all the pq part and reduced r<=s.
    else

       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(ntot,gmo,tmp,dimP,dimQ)&
       !$OMP PRIVATE(r,s,rs,ibatch,ipack,ncopy)
       do s=1,ntot
          do r=1,s
             rs = r + (s-1)*ntot
             ibatch = 1 + (rs-1)*dimP*dimQ
             ipack = 1 + (s*(s-1)/2 + r-1)*dimP*dimQ
             ncopy = dimP*dimQ
             call dcopy(ncopy,gmo(ibatch),1,tmp(ipack),1)
          end do
       end do
       !$OMP END PARALLEL DO

       ncopy = (ntot*(ntot-1)/2 + ntot)*dimP*dimQ
       ! accumulate tile
       if (nnod>1.and.pack_gmo%itype==TT_TILED_DIST) then
          call time_start_phase(PHASE_COMM)
          call tensor_accumulate_tile(pack_gmo,tile,tmp(1:ncopy),ncopy,lock_set=.true.)
          call time_start_phase(PHASE_WORK)
       else if (nnod>1.and.pack_gmo%itype==TT_TILED) then
          call daxpy(ncopy,1.0E0_realk,tmp,1,pack_gmo%ti(tile)%t(:),1)
       else
          call daxpy(ncopy,1.0E0_realk,tmp,1,pack_gmo%elm2(:,tile),1)
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
       else if (nnod>1.and.pack_gmo%itype==TT_TILED) then
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
       else if (nnod>1.and.pack_gmo%itype==TT_TILED) then
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

  subroutine get_mo_integral_par(integral,trafo1,trafo2,trafo3,trafo4,mylsitem,local,collective,order)
    implicit none
    type(tensor),intent(inout)   :: integral
    type(tensor),intent(inout)   :: trafo1,trafo2,trafo3,trafo4
    type(lsitem), intent(inout) :: mylsitem
    logical, intent(in) :: local
    logical, intent(inout) :: collective
    integer, intent(in), optional :: order(4)
    !Integral stuff
    integer :: alphaB,gammaB,dimAlpha,dimGamma
    integer :: dim1,dim2,dim3,MinAObatch
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
    integer :: iorb,nthreads,magic
    integer :: idx,nb,n1,n2,n3,n4,fa,fg,la,lg,i,k,myload,nba,nbg,biA,biG,bsA,bsG
#ifdef VAR_ICHOR
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
#else
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    Character(80)        :: FilenameCS,FilenamePS
    Character(80),pointer:: BatchfilenamesCS(:,:)
    Character(80),pointer:: BatchfilenamesPS(:,:)
    logical :: FoundInMem,doscreen
    integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
    TYPE(DECscreenITEM)  :: DecScreen
#endif
    integer, pointer :: batchdimAlpha(:),batchdimGamma(:)
    Character        :: INTSPEC(5)
    logical :: fullRHS
    integer :: MaxAllowedDimAlpha,MaxActualDimAlpha,nbatchesAlpha
    integer :: MaxAllowedDimGamma,MaxActualDimGamma,nbatchesGamma
    real(realk), pointer :: w1(:),w2(:)
    real(realk) :: MemFree,nrm
    integer(kind=long) :: maxsize
    logical :: master
    integer(kind=ls_mpik) :: me, nnod
    integer, pointer :: jobdist(:)
    real(realk), pointer :: work(:)
    integer(kind=long) :: w1size, w2size
    real(realk), parameter :: fraction_of = 0.8E0_realk

    call time_start_phase( PHASE_WORK )


    master  = .true.
    me      = 0
    nnod    = 1
    magic   = 3
#ifdef VAR_MPI
    master  = (infpar%lg_mynum == infpar%master)
    me      = infpar%lg_mynum
    nnod    = infpar%lg_nodtot
#endif


    nb = trafo1%dims(1)
    n1 = trafo1%dims(2)
    n2 = trafo2%dims(2)
    n3 = trafo3%dims(2)
    n4 = trafo4%dims(2)
    if( integral%dims(1) /= n1 .or. integral%dims(2) /= n2 .or. &
         & integral%dims(3) /= n3 .or. integral%dims(4) /= n4)then
       call lsquit("EEROR(get_mo_integral_par)wrong dimensions of the integrals&
            & or the transformation matrices",-1)
    endif

    ! Set integral info
    ! *****************
    INTSPEC(1)               = 'R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)               = 'R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)               = 'R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)               = 'R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)               = 'C' !C = Coulomb operator
#ifdef VAR_ICHOR
    iprint = 0           !print level for Ichor Integral code
    MoTrans = .FALSE.    !Do not transform to MO basis! 
    NoSymmetry = .FALSE. !Use Permutational Symmetry! 
    SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
    !Determine the full number of AO batches - not to be confused with the batches of AOs
    !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
    iAO = 1
    call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
    nullify(AOGammabatchinfo)
    nullify(AOalphabatchinfo)    
#else
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
#endif
    !==================================================
    !                  Batch construction             !
    !==================================================


    ! Get free memory and determine maximum batch sizes
    ! -------------------------------------------------
    if(master)then
#ifdef VAR_MPI
       call time_start_phase( PHASE_COMM )
       if(.not.local)call wake_slaves_for_simple_mo(integral,trafo1,trafo2,trafo3,&
            &trafo4,mylsitem,collective)
       call time_start_phase( PHASE_WORK )
#endif

#ifdef VAR_ICHOR
       !Determine the minimum allowed AObatch size MinAObatch
       !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
       !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
       !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
       !'R'  !Specifies that it is the Regular AO basis that should be batched
       iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
       call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
#else
       call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
#endif
       call get_currently_available_memory(MemFree)


       nba = nb
       nbg = nb
       alp: do i = MinAObatch, nb
          gamm: do k = MinAObatch, nb

             maxsize=max(max(nb**2*i*k,n1*n2*k*i),n1*n2*n3*n4)
             maxsize=maxsize + max(n1*nb*i*k,n1*n2*n3*k)
             if(collective) maxsize = maxsize + n1*n2*n3*n4

             if(float(maxsize*8)/(1024.0**3) > fraction_of*MemFree )then
                if(nba <= MinAObatch .and. nbg<= MinAObatch .and. collective)then
                   collective = .false.
                else
                   nba = i
                   nbg = k - 1
                   exit alp
                endif
             endif

          enddo gamm
       enddo alp


       if(DECinfo%manual_batchsizes)then
          nbg = max(DECinfo%ccsdGbatch,MinAObatch)
          nba = max(DECinfo%ccsdAbatch,MinAObatch)
       else
          if((nb/nba)*(nb/nbg)<magic*nnod.and.(nba>MinAObatch).and.nnod>1)then
             nba=(nb/(magic*nnod))
             if(nba<MinAObatch)nba=MinAObatch
          endif

          if((nb/nba)*(nb/nbg)<magic*nnod.and.(nba==MinAObatch).and.nnod>1)then
             do while((nb/nba)*(nb/nbg)<magic*nnod)
                nbg=nbg-1
                if(nbg<=MinAObatch)exit
             enddo
             if(nbg<MinAObatch)nbg=MinAObatch
          endif
       endif

       maxsize=max(max(nb**2*nba*nbg,n1*n2*nba*nbg),n1*n2*n3*n4)
       maxsize=maxsize + max(n1*nb*nba*nbg,n1*n2*n3*nbg)
       if(collective) maxsize = maxsize + n1*n2*n3*n4

       if(float(maxsize*8)/(1024.0**3) > fraction_of*MemFree)call lsquit("ERROR(get_mo_integral_par)not enough memory",-1)

       MaxAllowedDimGamma = nbg
       MaxAllowedDimAlpha = nba

    endif


    if(.not.local)then
       integral%access_type = AT_ALL_ACCESS
       trafo1%access_type   = AT_ALL_ACCESS
       trafo2%access_type   = AT_ALL_ACCESS
       trafo3%access_type   = AT_ALL_ACCESS
       trafo4%access_type   = AT_ALL_ACCESS
#ifdef VAR_MPI
       call time_start_phase( PHASE_COMM )
       call ls_mpibcast(MaxAllowedDimAlpha,infpar%master,infpar%lg_comm)
       call ls_mpibcast(MaxAllowedDimGamma,infpar%master,infpar%lg_comm)
       call time_start_phase( PHASE_WORK )
#endif
    endif

    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
#ifdef VAR_ICHOR
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
#else
    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchGamma,nb)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
         & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
         &nbatchesGamma,orb2BatchGamma,'R')
#endif

    if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma,&
         & 'with maximum size',MaxActualDimGamma

#ifndef VAR_ICHOR
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
#endif

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

#ifdef VAR_ICHOR
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
#else
    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nb)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
         & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')
#endif

    if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha&
         &, 'with maximum size',MaxActualDimAlpha

#ifndef VAR_ICHOR
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
#endif

    maxsize=max(max(nb**2*MaxActualDimAlpha*MaxActualDimGamma,n1*n2*MaxActualDimAlpha*MaxActualDimGamma),n1*n2*n3*n4)
    w1size = maxsize
    call mem_alloc( w1, w1size )
    maxsize=max(n1*nb*MaxActualDimAlpha*MaxActualDimGamma,n1*n2*n3*MaxActualDimGamma)
    w2size = maxsize
    call mem_alloc( w2, w2size )
    if(collective)then
       call mem_alloc(work,(i8*n1)*n2*n3*n4)
       work = 0.0E0_realk
    endif


    ! ************************************************
    ! *  precalculate the full schreening matrix     *
    ! ************************************************
#ifdef VAR_ICHOR
     !Calculate Screening integrals 
     SameMOL = .TRUE. !Specifies same molecule on all centers 
     call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
#else
    ! This subroutine builds the full screening matrix.
    call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
         & nbatchesAlpha,nbatchesGamma,INTSPEC)
    IF(mylsitem%setting%scheme%cs_screen .OR. mylsitem%setting%scheme%ps_screen)THEN
       call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
            & nb,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
       call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
            & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
            & batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
    ENDIF
#endif

    call mem_alloc(jobdist,nbatchesAlpha*nbatchesGamma)
    !JOB distribution
#ifdef VAR_MPI


#ifdef VAR_ICHOR
    call mem_alloc(batchdimAlpha,nbatchesAlpha)
    do idx=1,nbatchesAlpha
       batchdimAlpha(idx) = AOAlphabatchinfo(idx)%dim 
    enddo
    call mem_alloc(batchdimGamma,nbatchesGamma)
    do idx=1,nbatchesGamma
       batchdimGamma(idx) = AOGammabatchinfo(idx)%dim 
    enddo
#endif
    call distribute_mpi_jobs(jobdist,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
         &batchdimGamma,myload,nnod,me)
#ifdef VAR_ICHOR
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchdimGamma)
#endif
#else
    jobdist = 0
#endif

    !print *,me,"has",batchindexGamma,batchindexAlpha
    !call lsmpi_barrier(infpar%lg_comm)

    myload = 0
    fullRHS = nbatchesGamma.EQ.1.AND.nbatchesAlpha.EQ.1

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches

#ifdef VAR_ICHOR
       lg = AOGammabatchinfo(gammaB)%dim               ! Dimension of gamma batch
       fg = AOGammabatchinfo(gammaB)%orbstart          ! First orbital index in gamma batch
       GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
       AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
       AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
#else
       lg  = batchdimGamma(gammaB)                     ! Dimension of gamma batch
       fg  = batch2orbGamma(gammaB)%orbindex(1)        ! First index in gamma batch
       biG = batchindexGamma(gammaB)
       bsG = batchsizeGamma(gammaB)
#endif
       BatchAlpha: do alphaB = 1, nbatchesAlpha

#ifdef VAR_ICHOR
          la = AOAlphabatchinfo(alphaB)%dim               ! Dimension of alpha batch
          fa = AOAlphabatchinfo(alphaB)%orbstart          ! First orbital index in alpha batch
          AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
          AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
          AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
#else
          la  = batchdimAlpha(alphaB)                              ! Dimension of alpha batch
          fa  = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
          biA = batchindexAlpha(alphaB)
          bsA = batchsizeAlpha(alphaB)
#endif
          !print '(I3,"have",8I7)',me,lg,fg,biG,bsG,la,fa,biA,bsA
          !call lsmpi_barrier(infpar%lg_comm)

          if( me /= jobdist(gammaB + (alphaB-1) *nbatchesGamma) ) cycle BatchAlpha

          if(DECinfo%PL>2)write (*, '("Rank",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")')&
               &me,alphaB,nbatchesAlpha,gammaB,nbatchesGamma

          myload     = myload + la * lg

#ifdef VAR_ICHOR
          call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nb,nb,la,lg,&
               & w1,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
               & AOGammaStart,AOGammaEnd,MoTrans,nb,nb,la,lg,NoSymmetry)
#else
          IF(doscreen) Mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
          IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p

          call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, Mylsitem%setting, w1,biA,&
               &biG,bsA,bsG,nb,nb,la,lg,fullRHS,INTSPEC)
#endif

          !something more sophisticated can be implemented here
          call dgemm('t','n',nb*la*lg,n1,nb,1.0E0_realk,w1,nb,trafo1%elm1,nb,0.0E0_realk,w2,nb*la*lg)
          call dgemm('t','n',la*lg*n1,n2,nb,1.0E0_realk,w2,nb,trafo2%elm1,nb,0.0E0_realk,w1,la*lg*n1)
          call dgemm('t','n',lg*n1*n2,n3,la,1.0E0_realk,w1,la,trafo3%elm1(fa),nb,0.0E0_realk,w2,lg*n1*n2)

          if(collective) then
             call dgemm('t','n',n1*n2*n3,n4,lg,1.0E0_realk,w2,lg,trafo4%elm1(fg),nb,1.0E0_realk,work,n1*n2*n3)
          else
             call dgemm('t','n',n1*n2*n3,n4,lg,1.0E0_realk,w2,lg,trafo4%elm1(fg),nb,0.0E0_realk,w1,n1*n2*n3)

             call time_start_phase( PHASE_COMM )
             call tensor_add(integral,1.0E0_realk,w1,wrk=w2,iwrk=maxsize, order = order)
             call time_start_phase( PHASE_WORK )
          endif

       enddo BatchAlpha
    enddo BatchGamma

    ! Free integral stuff
    ! *******************
#ifdef VAR_ICHOR
    call FREE_SCREEN_ICHORERI()
    call mem_dealloc(AOGammabatchinfo)
    call mem_dealloc(AOAlphabatchinfo)
#else
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
#endif

    call mem_dealloc(jobdist)


    call mem_dealloc( w1 )
    call mem_dealloc( w2 )

#ifdef VAR_MPI
    call time_start_phase( PHASE_IDLE )
    call lsmpi_barrier(infpar%lg_comm)
    if(collective)then
       call time_start_phase( PHASE_COMM )
       call lsmpi_allreduce(work,(i8*n1)*n2*n3*n4,infpar%lg_comm)
       call tensor_convert(work,integral, order = order )
    endif
    call time_start_phase( PHASE_WORK )
#else
    call tensor_convert(work,integral, order = order )
#endif
    if(collective) call mem_dealloc( work )

    if(DECinfo%PL>2)then
       call print_norm(integral,nrm)
       if(master) print *," NORM of the integral :",nrm
    endif

    if(.not.local)then
       integral%access_type = AT_MASTER_ACCESS
       trafo1%access_type = AT_MASTER_ACCESS
       trafo2%access_type = AT_MASTER_ACCESS
       trafo3%access_type = AT_MASTER_ACCESS
       trafo4%access_type = AT_MASTER_ACCESS
    endif


  end subroutine get_mo_integral_par

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
  !> CC model:
  integer ::  ccmodel
  !> SCF transformation matrices:
  real(realk), pointer  :: Co(:,:), Cv(:,:)
  !> performed MO-based CCSD calculation ?
  logical :: mo_ccsd
  !> array with gmo on output:
  type(tensor) :: pgmo_diag, pgmo_up, govov
  !> variables used for MO batch and integral transformation
  type(MObatchInfo) :: MOinfo
  !> LS item information
  type(lsitem) :: MyLsItem


  call mpi_communicate_get_gmo_data(mo_ccsd,MyLsItem,Co,Cv, &
       & pgmo_diag,pgmo_up,nb,no,nv,nbatch,ccmodel)

  ! the slave call the routine to get MO int.
  call get_t1_free_gmo(mo_ccsd,MyLsItem,Co,Cv,govov, &
       & pgmo_diag,pgmo_up,nb,no,nv,ccmodel,MOinfo)

  ! deallocate slave stuff:
  call mem_dealloc(Co)
  call mem_dealloc(Cv)
  call ls_free(MyLsItem)
  if (ccmodel/=MODEL_RPA) then 
     call mem_dealloc(MOinfo%dimInd1)
     call mem_dealloc(MOinfo%dimInd2)
     call mem_dealloc(MOinfo%StartInd1)
     call mem_dealloc(MOinfo%StartInd2)
     call mem_dealloc(MOinfo%dimTot)
     call mem_dealloc(MOinfo%tileInd)
  end if

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

  call wake_slaves_for_simple_mo(integral,trafo1,trafo2,trafo3,trafo4,mylsitem,c)
  call get_mo_integral_par(integral,trafo1,trafo2,trafo3,trafo4,mylsitem,.false.,c)
  call ls_free(mylsitem)

end subroutine get_mo_integral_par_slave
#endif
