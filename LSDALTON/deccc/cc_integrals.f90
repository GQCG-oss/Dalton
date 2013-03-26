!> @file
!> Simple integrals related
!> \author Marcin Ziolkowski
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
  use Integralparameters
  use integralinterfaceMOD

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use array2_simple_operations
  use array4_simple_operations

  interface getL
     module procedure getL_simple
     module procedure getL_simple_from_gmo
     module procedure getL_diff
  end interface

contains

  !> \brief Get full two-electron integrals in AO basis using DALTONs integral program
  subroutine get_full_eri(mylsitem,nbasis,g_ao)

    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer, intent(in) :: nbasis
    type(array4),intent(inout) :: g_ao
    type(matrix) :: g_matrix
    integer, dimension(4) :: ao_dims

    ao_dims=[nbasis,nbasis,nbasis,nbasis]

    write(DECinfo%output,'(a)') 'info :: calculating two-electron integrals'
    g_ao = array4_init_standard(ao_dims)
    ! KK Quick fix: Filename associated with g_ao
    g_ao%filename = 'gao'

    call ii_get_4center_eri(DECinfo%output,DECinfo%output,&
         & mylsitem%setting,g_ao%val,nbasis,nbasis,nbasis,nbasis)

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
    if(DECinfo%show_time) write(DECinfo%output,'(a,f16.3,a)') &
         ' time :: integral transformation : ',endtime-starttime,' s'

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

  !> \brief Simple routine to transform integrals to (AI|BJ) for storing all AOs in memory
  !> \author Marcin Ziolkowski
  !> \param gao Full matrix with AO integrals
  !> \param left1 First index transformation (to A)
  !> \param left2 Second index transformation (to I)
  !> \param right1 Third index transformation (to B)
  !> \param right2 Fourth index transformation (to J)
  !> \return Array4 with intrgrals (AI|BJ)
  function get_gmo_vovo(gao,left1,left2,right1,right2) result(gmo)

    implicit none
    type(array4), intent(inout) :: gao
    type(array4) :: gmo
    type(array4) :: tmp1,tmp2,tmp3
    type(array2), intent(inout) :: left1,left2,right1,right2
    real(realk), pointer :: AA(:,:), BB(:,:), CC(:,:), DD(:,:)
    real(realk) :: starttime,endtime
    integer, dimension(4) :: dims,rotate
    integer :: nbas,nocc,nvirt
    integer :: a,b,c,d,p,q,r,s

    call cpu_time(starttime)

    nocc  = left2%dims(2)
    nvirt = left1%dims(2)
    nbas = left1%dims(1)

    rotate=[2,3,4,1]

    ! transform first index to occupied
    call array4_read(gao)
    tmp1 = array4_init([nocc,nbas,nbas,nbas])
    call array4_contract1(gao,left2,tmp1,.true.) ! gao[MuNu,AlBe] -> tmp1[JNu,AlBe]

    ! deallocate ao integrals
    call array4_dealloc(gao)

    call array4_reorder(tmp1,[3,4,2,1]) ! tmp1[JNu,AlBe] -> tmp1[AlBe,NuJ]
    tmp2 = array4_init([nocc,nbas,nbas,nocc])
    call array4_contract1(tmp1,right2,tmp2,.true.) ! tmp1[AlBe,NuJ] -> tmp2[IBe,NuJ]
    print *,'half trans :',tmp2*tmp2
    call array4_free(tmp1)

    call array4_reorder(tmp2,[2,1,3,4]) ! tmp2[IBe,NuJ] -> tmp2[BeI,NuJ]
    tmp3 = array4_init([nvirt,nocc,nbas,nocc])
    call array4_contract1(tmp2,right1,tmp3,.true.) ! tmp2[BeI,NuJ] -> tmp3[AI,NuJ]
    call array4_free(tmp2)
    print *,'AI,NuJ ',tmp3*tmp3

    call array4_reorder(tmp3,[3,4,1,2]) ! tmp3[AI,NuJ] -> tmp3[NuJ,AI]
    gmo = array4_init([nvirt,nocc,nvirt,nocc])
    call array4_contract1(tmp3,left1,gmo,.true.) ! tmp3[NuJ,AI] -> gmo[BJ,AI]
    call array4_free(tmp3)

    call cpu_time(endtime)
    if(DECinfo%show_time) write(DECinfo%output,'(a,f16.3,a)') &
         ' time :: integral transformation : ',endtime-starttime,' s'

    return
  end function get_gmo_vovo


  !> \brief Carry out 2-electron Fock transformation on matrix U, i.e.
  !> FockU = 2J(U) - K(U)
  !> where J(U) and K(U) are the coulomb and exchange transformations, respectively.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_fock_transformation(FockU,U,MyLsitem,symmetric)

    implicit none
    !> Fock transformation on U matrix
    type(matrix), intent(inout) :: FockU
    !> Matrix to carry out Fock transformation on
    type(matrix), intent(in) :: U
    !> LS DALTON info
    type(lsitem), intent(inout) :: MyLsItem
    !> Is U symmetric (true) or not (false)?
    logical, intent(in) :: symmetric
    real(realk) :: tcpu1,twall1,tcpu2,twall2,tcpu,twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! Sanity check
    if(U%nrow /= U%ncol) then
       call lsquit('dec_fock_transformation:&
            & Matrix U must be quadratic!',-1)
    end if

    ! Carry out Fock transformation on U
    call II_get_Fock_mat(DECinfo%output, DECinfo%output, &
         & MyLsitem%setting,U,symmetric,FockU,1,.FALSE.)

    call LSTIMER('DEC: FOCK TRANS',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    DECinfo%integral_time_cpu = DECinfo%integral_time_cpu + (tcpu2-tcpu1)
    DECinfo%integral_time_wall = DECinfo%integral_time_wall + (twall2-twall1)

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
    integer :: dummy
    TYPE(DECscreenITEM)   :: DecScreen
    logical :: doscreen
    integer :: ndim(4),i
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

    ! Set integral screening
    dummy=0
    doscreen = mylsitem%setting%scheme%cs_screen .OR. mylsitem%setting%scheme%ps_screen
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mylsitem%setting,dummy,1,1,intspecConvert)
    IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabLHS
    IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabRHS

    ! Get AO integrals
    call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output,Mylsitem%SETTING,&
         & gao,1,1,ndim(3),ndim(4),ndim(1),ndim(2),ndim(3),ndim(4),.true.,1,intSpecConvert)
    call free_decscreen(DECSCREEN)
    nullify(mylsitem%setting%LST_GAB_RHS)
    nullify(mylsitem%setting%LST_GAB_LHS)

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


end module ccintegrals
