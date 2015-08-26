!> @file
!> Routines for matrix type and operations compatibile with array4 type
!> \author Marcin Ziolkowski
!> \todo Document this file better
module array2_simple_operations

  ! Outside DEC directory
  use precision
  use dec_typedef_module
  use memory_handling!, only: mem_alloc, mem_dealloc


  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils
  use dec_tools_module

  integer :: CreatedArrays2=0
  integer :: DestroyedArrays2=0



  !> Overloaded constructor
  interface array2_init
    module procedure array2_init_plain
    module procedure array2_init_data
  end interface array2_init

  !> Overloaded function for adding arrays
  interface array2_add
    module procedure array2_add_plain
    module procedure array2_add_scale
  end interface array2_add

  !> Overloaded operator for addition
  interface operator(+)
    module procedure array2_add_plain
  end interface

  !> Overloaded operator for dot-product
  interface operator(*)
    module procedure array2_dotproduct
  end interface

  contains

  !> \brief Initialize empty array
  function array2_init_plain(dims) result(array)

    implicit none
    type(array2) :: array
    integer, dimension(2), intent(in) :: dims

    array%dims=dims
    if(associated(array%val)) nullify(array%val)
    call mem_alloc(array%val,dims(1),dims(2) )
    array%val = 0.0E0_realk

    CreatedArrays2=CreatedArrays2+1

    return
  end function array2_init_plain

  !> \brief Initialize an array from fortran matrix
  function array2_init_data(dims,matrix) result(array)

    implicit none
    type(array2) :: array
    integer, dimension(2), intent(in) :: dims
    real(realk), dimension(:,:), intent(in) :: matrix

    array%dims=dims
    if(associated(array%val)) nullify(array%val)
    call mem_alloc(array%val,dims(1),dims(2) )

    CreatedArrays2=CreatedArrays2+1

    ! Assign
    array%val=matrix

    return
  end function array2_init_data

  !> \brief Get a matrix by similarity transformation
  function array2_similarity_transformation(left,matrix,right,final_dims) &
      result(final)

    implicit none
    type(array2), intent(inout) :: left,right,matrix
    integer, dimension(2), intent(in) :: final_dims
    type(array2) :: final, tmp


    tmp = array2_init([left%dims(2),matrix%dims(2)])
    final=array2_init(final_dims)

    call array2_matmul(left,matrix,tmp,'t','n',1.0E0_realk,0.0E0_realk)
    call array2_matmul(tmp,right,final,'n','n',1.0E0_realk,0.0E0_realk)

    call array2_free(tmp)

    return
  end function array2_similarity_transformation

  ! KK.
  ! get a matrix by similarity transformation,
  ! from MO basis to AO basis:
  ! final = left * matrix * right^T
  function array2_similarity_transformation_MO_to_AO(left,matrix,right,final_dims) &
      result(final)

    implicit none
    type(array2), intent(inout) :: left,right,matrix
    integer, dimension(2), intent(in) :: final_dims
    type(array2) :: final, tmp


    tmp = array2_init([left%dims(1),matrix%dims(2)])
    final=array2_init(final_dims)

    call array2_matmul(left,matrix,tmp,'n','n',1.0E0_realk,0.0E0_realk)
    call array2_matmul(tmp,right,final,'n','t',1.0E0_realk,0.0E0_realk)

    call array2_free(tmp)

    return
  end function array2_similarity_transformation_MO_to_AO


  !> \brief Duplicate array
  function array2_duplicate(this) result(newarray)

    implicit none
    type(array2), intent(in) :: this
    type(array2) :: newarray

    newarray = array2_init(this%dims,this%val)

  end function array2_duplicate

  !> \brief Destroy array
  subroutine array2_free(array)

    implicit none
    type(array2), intent(inout) :: array

    if(associated(array%val)) then
      call mem_dealloc(array%val)
      array%val => null()
    end if

    DestroyedArrays2 = DestroyedArrays2 + 1

    return
  end subroutine array2_free

  !> \brief Print statistics of this simple array routines
  subroutine array2_print_statistics(output)

    implicit none
    integer, intent(in) :: output

    write(DECinfo%output,'(/,a)')    ' -------------------------------'
    write(DECinfo%output,'(a)')      '        Array2 statistics       '
    write(DECinfo%output,'(a)')      ' -------------------------------'
    write(DECinfo%output,'(a,i4)')   'Number of created arrays   : ',CreatedArrays2
    write(DECinfo%output,'(a,i4)')   'Number of destroyes arrays : ',DestroyedArrays2
    write(DECinfo%output,'(a,i4,/)') 'Orphaned arrays            : ',CreatedArrays2-DestroyedArrays2

    return
  end subroutine array2_print_statistics

  !> \brief Transpose array
  subroutine array2_transpose(array)

    implicit none
    type(array2), intent(inout) :: array
    real(realk), pointer :: dat(:,:) => null()
    integer, dimension(2) :: dims
    integer :: i,j

    dims(1) = array%dims(2)
    dims(2) = array%dims(1)

    ! transpose
    call mem_alloc(dat,dims(1),dims(2) )
    do i=1,dims(1)
      do j=1,dims(2)
        dat(i,j) = array%val(j,i)
      end do
    end do
    call mem_dealloc(array%val)

    ! assign transposed data
    array%val => dat
    array%dims(1) = dims(1)
    array%dims(2) = dims(2)

    return
  end subroutine array2_transpose

  !> \brief Take symmetric component of Asym = 1/2*(A + A^T)
  function array2_symmetric_component(A) result(Asym)

    implicit none
    type(array2), intent(in) :: A
    type(array2) :: Asym
    type(array2) :: At
    real(realk), pointer :: dat(:,:) => null()
    integer :: dim1, dim2
    integer :: i,j

    dim1 = A%dims(2)
    dim2 = A%dims(1)
    if(dim1 /= dim2) then
       call lsquit('array2_symmetric_component: dim1 /= dim2',DECinfo%output)
    end if

    At=array2_init([dim1,dim2],A%val)
    call array2_transpose(At)

    Asym = array2_add(A,At)
    Asym%val = 0.5E0_realk*Asym%val

    call array2_free(At)

  end function array2_symmetric_component


  !> \brief Print some info
  subroutine array2_info(array,output)

    implicit none
    type(array2), intent(in) :: array
    integer, intent(in) :: output
    real(realk) :: memory

    write(DECinfo%output,'(/,a)')   '    Array2 Info     '
    write(DECinfo%output,'(a)')     '--------------------'
    write(DECinfo%output,'(a,i4)')  '    Dim1 : ',array%dims(1)
    write(DECinfo%output,'(a,i4)')  '    Dim2 : ',array%dims(2)

    memory = dble(array%dims(1))*dble(array%dims(2))*realk/dble(2**20)
    write(DECinfo%output,'(a,f10.2,a)') '  Memory : ',memory,' MB'

    return
  end subroutine array2_info

  !> \brief Dot product of two matrices
  function array2_dotproduct(a,b) result(res)

    implicit none
    type(array2), intent(in) :: a,b
    real(realk) :: ddot
    real(realk) :: res
    integer :: nelements
    integer :: i

    res=0.0E0_realk

!#ifdef USE_BLAS
    nelements = a%dims(1)*a%dims(2)
    res = ddot(nelements,a%val,1,b%val,1)
!#else
!    do i=1,a%dims(2)
!      res = res + dot_product(a%val(:,i),b%val(:,i))
!    end do
!#endif

    return
  end function array2_dotproduct

  !> \brief Get a square 2norm (aka 1norm from cc_driver)
  function array2_norm(array) result(val)

    implicit none
    type(array2), intent(in) :: array
    real(realk) :: val

    val = 0.0E0_realk
    val = array2_dotproduct(array,array)
    return
  end function array2_norm

  !> \brief Set array to a given value
  subroutine array2_set_value(array,val)

    implicit none
    type(array2), intent(inout) :: array
    real(realk), intent(in) :: val

    if(associated(array%val)) then
      array%val=val
    end if

    return
  end subroutine array2_set_value

  !> \brief Set array to zero
  subroutine array2_zero(array)

    implicit none
    type(array2), intent(inout) :: array
    call array2_set_value(array,0.0E0_realk)
    return
  end subroutine array2_zero

  !> \brief Add two arrays
  function array2_add_plain(a,b) result(c)

    implicit none
    type(array2), intent(in) :: a,b
    type(array2) :: c
    integer, dimension(2) :: dims

    dims(1) = a%dims(1)
    dims(2) = a%dims(2)
    c = array2_init(dims)
    c%val = a%val + b%val
    return
  end function array2_add_plain

  !> \brief Add two arrays with scaling them
  function array2_add_scale(alpha,a,beta,b) result(c)

    implicit none
    real(realk), intent(in) :: alpha,beta
    type(array2), intent(in) :: a,b
    type(array2) :: c
    integer, dimension(2) :: dims

    dims(1) = a%dims(1)
    dims(2) = a%dims(2)
    c = array2_init(dims)
    c%val = alpha*a%val + beta*b%val
    return
  end function array2_add_scale

  !> \brief Add to one matrix (daxpy) a+=beta*b
  subroutine array2_add_to(a,beta,b)

    implicit none
    type(array2), intent(in) :: b
    type(array2) :: a
    real(realk), intent(in) :: beta

    a%val=a%val+beta*b%val
    return
  end subroutine array2_add_to

  !> \brief Copy b to a
  subroutine array2_copy(a,b)

    implicit none
    type(array2), intent(inout) :: a
    type(array2), intent(in) :: b

    if(a%dims(1) /= b%dims(1) .or. a%dims(2) /= b%dims(2) ) &
      stop 'error :: array2_copy -> dimensions dont match'

    a%val = b%val

    return
  end subroutine array2_copy

  !> \brief Matmul -> c=beta*c+alpha*matmul(op(a),op(b))
  subroutine array2_matmul(a,b,c,trans_a,trans_b,alpha,beta)

    implicit none
    type(array2), intent(inout) :: a,b,c
    character, intent(in) :: trans_a, trans_b
    real(realk), intent(in) :: alpha, beta

    ! transpose
    if(trans_a == 't' .or. trans_a == 'T') call array2_transpose(a)
    if(trans_b == 't' .or. trans_b == 'T') call array2_transpose(b)

!#ifdef USE_BLAS
    call dgemm('n','n',a%dims(1),b%dims(2),a%dims(2),alpha,a%val, &
            a%dims(1),b%val,b%dims(1),beta,c%val,c%dims(1))
!!$#else
!!$    c%val = beta*c%val+alpha*matmul(a%val,b%val)
!!$#endif

    ! transpose back
    if(trans_a == 't' .or. trans_a == 'T') call array2_transpose(a)
    if(trans_b == 't' .or. trans_b == 'T') call array2_transpose(b)

    return
  end subroutine array2_matmul

  !> \brief Print
  subroutine array2_print(this,label)

    implicit none
    type(array2), intent(in) :: this
    character(*), intent(in) :: label
    integer :: i,j

    write(DECinfo%output,*) 'ARRAY2 LABEL: ', label
    call array2_info(this,DECinfo%output)

    write(DECinfo%output,'(a)') ' -- Data --'
    write(DECinfo%output,'(a)') ' -> i,j,data(i,j)'
    do i=1,this%dims(1)
      do j=1,this%dims(2)
        write(DECinfo%output,'(2i4,f16.10)') i,j,this%val(i,j)
      end do
    end do

    return
  end subroutine array2_print

  !> \brief Unit matrix
  function array2_unit(dims) result(new)

    implicit none
    type(array2) :: new
    integer, dimension(2), intent(in) :: dims
    integer :: i

    if( dims(1) /= dims(2) ) &
          stop ' error :: cannot make not-square unit matrix'

    new = array2_init(dims)
    do i=1,dims(1)
      new%val(i,i) = 1.0E0_realk
    end do

    return
  end function array2_unit



  !> \brief Extract EOS indices from two-dimensional array.
  !> It is possible to extract EOS indices from either columns or rows, and
  !> and using either the occupied or the virtual partitioning scheme.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine array2_extract_EOS(AOS,MyFragment,occ_or_virt,row_or_column,EOS)

    implicit none
    !> Matrix where rows, columns or both are AOS indices
    type(array2),intent(in) :: AOS
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    !> Extract occupied ('o' or 'O') or virtual ('v' or 'V') EOS indices from AOS matrix
    character(len=1), intent(in) :: occ_or_virt
    !> Extract EOS indices from rows ('r' or 'R') or columns ('c' or 'C') of AOS matrix
    character(len=1), intent(in) :: row_or_column
    !> Output matrix where EOS indices have been extracted (also initialized here)
    type(array2),intent(inout) :: EOS
    integer :: dim1,dim2
    integer :: i, j,idx1, idx2, ES
    logical :: row,occ, something_wrong



    ! Check input
    ! ***********

    ! Extract occupied or virtual EOS indices?
    if( (occ_or_virt=='o') .or. (occ_or_virt=='O') ) then
       occ=.true.
    elseif( (occ_or_virt=='v') .or. (occ_or_virt=='V') ) then
       occ=.false.
    else
       call lsquit('array2_extract_EOS: occ_or_virt argument must be O or V!',DECinfo%output)
    end if

    ! Extract rows or columns?
    if( (row_or_column=='r') .or. (row_or_column=='R') ) then
       row=.true.
    elseif( (row_or_column=='c') .or. (row_or_column=='C') ) then
       row=.false.
    else
       call lsquit('array2_extract_EOS: row_or_column argument must be R or C!',DECinfo%output)
    end if


    ! Sanity checks and set extraction scheme
    ! ***************************************
    ! ES=1: Extract occupied row EOS
    ! ES=2: Extract occupied column EOS
    ! ES=3: Extract virtual row EOS
    ! ES=4: Extract virtual column EOS

    something_wrong=.false.
    if(occ .and. row ) then
       ! Number of rows in input matrix must correspond to occupied AOS
       if(AOS%dims(1) /= MyFragment%noccAOS) something_wrong=.true.
       ES=1  ! extraction scheme 1
    end if
    if(occ .and. (.not. row) ) then
       ! Number of columns in input matrix must correspond to occupied AOS
       if(AOS%dims(2) /= MyFragment%noccAOS) something_wrong=.true.
       ES=2  ! extraction scheme 2
    end if
    if( (.not. occ) .and. row ) then
       ! Number of rows in input matrix must correspond to virtual AOS
       if(AOS%dims(1) /= MyFragment%nvirtAOS) something_wrong=.true.
       ES=3  ! extraction scheme 3
    end if
    if( (.not. occ) .and. (.not. row) ) then
       ! Number of columns in input matrix must correspond to virtual AOS
       if(AOS%dims(2) /= MyFragment%nvirtAOS) something_wrong=.true.
       ES=4  ! extraction scheme 4
    end if
    if(something_wrong) then
       write(DECinfo%output,*) 'Input matrix dims = ', AOS%dims
       write(DECinfo%output,*) 'Occupied/virtual AOS = ', MyFragment%noccAOS, &
            & MyFragment%nvirtAOS
       write(DECinfo%output,*) 'Logical row/occupied', row,occ
       call lsquit('Inconsistent input in array2_extract_EOS',DECinfo%output)
    end if



    ! Initialize output EOS array
    ! ***************************
    select case(ES)
    case(1)  ! ES=1: Extract occupied row EOS, same columns as AOS input
       dim1 = MyFragment%noccEOS
       dim2 = AOS%dims(2)
    case(2)  ! ES=2: Extract occupied column EOS, same rows as AOS input
       dim1 = AOS%dims(1)
       dim2 = MyFragment%noccEOS
    case(3)  ! ES=3: Extract virtual row EOS, same columns as AOS input
       dim1 = MyFragment%nvirtEOS
       dim2 = AOS%dims(2)
    case(4)  ! ES=4: Extract virtual column EOS, same rows as AOS input
       dim1 = AOS%dims(1)
       dim2 = MyFragment%nvirtEOS
    end select
    EOS = array2_init([dim1,dim2])


    ! Loop over orbitals and extract EOS orbitals
    ! *******************************************
    do j=1,dim2
       do i=1,dim1

          ! Extraction scheme
          ! -----------------
          select case(ES)
          case(1)  ! ES=1: Extract occupied row EOS
             idx1=MyFragment%idxo(i)
             idx2=j
          case(2)  ! ES=2: Extract occupied column EOS
             idx1=i
             idx2=MyFragment%idxo(j)
          case(3)  ! ES=3: Extract virtual row EOS
             idx1=MyFragment%idxu(i)
             idx2=j
          case(4)  ! ES=4: Extract virtual column EOS
             idx1=i
             idx2=MyFragment%idxu(j)
          end select

          ! Set EOS matrix
          EOS%val(i,j) = AOS%val(idx1,idx2)

       end do
    end do


  end subroutine array2_extract_EOS




  !> \brief Create MO coefficient matrix with occupied indices for only the EOS space,
  !> i.e. only occupied MOs assigned to the central atom in the fragment
  !> - or to one of the two atoms in case of a pair fragment.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine extract_occupied_EOS_MO_indices(Ceos,MyFragment)

    implicit none
    !> MO coeffiecient matrix for fragment occupied EOS
    type(array2),intent(inout) :: Ceos
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    ! Number of basis functions in fragment
    integer :: nbasis
    ! Number of occupied orbitals in EOS
    integer :: nocc_eos
    integer :: i, idx


    nocc_eos = MyFragment%noccEOS
    nbasis = MyFragment%nbasis

    ! Loop over orbitals and extract EOS orbitals
    do i=1,nocc_eos
       idx=MyFragment%idxo(i)
       Ceos%val(1:nbasis,i) = MyFragment%Co(1:nbasis,idx)
    end do


  end subroutine extract_occupied_EOS_MO_indices



  !> \brief Create MO coefficient matrix with virtual indices for only the EOS space,
  !> i.e. only virtual MOs assigned to the central atom in the fragment
  !> - or to one of the two atoms in case of a pair fragment.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine extract_virtual_EOS_MO_indices(Ceos,MyFragment)

    implicit none
    !> MO coeffiecient matrix for fragment virtual EOS
    type(array2),intent(inout) :: Ceos
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    ! Number of basis functions in fragment
    integer :: nbasis
    ! Number of occupied orbitals in EOS
    integer :: nvirt_eos
    integer :: a, ax

    nvirt_eos = MyFragment%nvirtEOS
    nbasis = MyFragment%nbasis

    ! Loop over orbitals and extract virtual EOS orbitals
    do a=1,nvirt_eos
       ax=MyFragment%idxu(a)
       Ceos%val(1:nbasis,a) = MyFragment%Cv(1:nbasis,ax)
    end do


  end subroutine extract_virtual_EOS_MO_indices



  !> \brief Construct transformation matrices used in MP2 integral scheme:
  !> 1. Transformation matrix between diagonal basis (AOS) and localized basis (AOS).
  !> 2. Transformation matrix from atomic basis (atomic extent) to Fock-diagonal basis (AOS).
  !> Matrices are constructed separately for occupied and virtual spaces,
  !> and all dimensions refer to a fragment (not the full molecule).
  !> Note: The array2 structures are also initialized here according to the dimensions in MyFragment.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine get_MP2_integral_transformation_matrices(MyFragment,CDIAGocc, CDIAGvirt, Uocc, Uvirt, &
       & EVocc, EVvirt)

    implicit none
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    !> Transforming from AO to occupied orbitals in the basis where the fragment Fock matrix is diagonal
    type(array2),intent(inout) :: CDIAGocc
    !> Transforming from AO to virtual orbitals in the basis where the fragment Fock matrix is diagonal
    type(array2),intent(inout) :: CDIAGvirt
    !> Occupied space: Transforming between local basis and basis where the fragment Fock matrix is diagonal
    !> - indices: (local,semi-canonical)
    type(array2),intent(inout) :: Uocc
    !> Virtual space: Transforming between local basis and basis where the fragment Fock matrix is diagonal
    !> - indices: (local,semi-canonical)
    type(array2),intent(inout) :: Uvirt
    !> Eigenvalues from diagonalization of occ-occ Fock matrix block
    real(realk), dimension(MyFragment%noccAOS),intent(inout) :: EVocc
    !> Eigenvalues from diagonalization of virt-virt Fock matrix block
    real(realk), dimension(MyFragment%nvirtAOS),intent(inout) :: EVvirt
    real(realk), pointer :: S(:,:)
    integer :: nocc,nvirt,nbasis,i
    integer, dimension(2) :: occocc,virtvirt,occAO,virtAO
    type(array2) :: Ctmp


    ! Init stuff
    ! **********
    nocc = MyFragment%noccAOS      ! (occ AOS size)
    nvirt = MyFragment%nvirtAOS   ! (virt AOS size)
    nbasis = MyFragment%nbasis        ! (atomic extent)
    occocc = [nocc,nocc]
    virtvirt = [nvirt,nvirt]
    occAO = [nbasis,nocc]
    virtAO = [nbasis,nvirt]



    ! ****************************************************************************************
    ! *                                  OCCUPIED ORBITAL SPACE                              *
    ! ****************************************************************************************


    ! 1. Solve Fock eigenvalue problem
    ! ********************************

    ! Eigenvectors
    Uocc = array2_init(occocc)

    ! The overlap matrix is simply the unit matrix because the local basis functions are orthogonal.
    call mem_alloc(S,nocc,nocc)
    S=0.0E0_realk
    do i=1,nocc
       S(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem

    call solve_eigenvalue_problem(nocc,MyFragment%ppfock,S,EVocc,Uocc%val)

    call mem_dealloc(S)

    ! Now Uocc contains the transformation coefficients between the local and diagonal basis.
    ! In particular, using capital/small letters for diagonal/local basis we have:
    ! phi_{K} = sum_{k} phi_{k} U_{k K}
    ! phi_{k} = sum_{K} phi_{K} U_{k K}



    ! 2. From AO to diagonal basis
    ! ****************************

    ! The transformation from the AO basis to the diagonal basis CDIAGocc is then:
    !
    ! CDIAGocc = CLOCALocc Uocc     (*)
    !
    ! where CLOCALocc are the basis MO coefficients in MyFragment%Co

    ! Local MO coefficients
    Ctmp = array2_init(occAO,MyFragment%Co)

    ! Determine CDIAGocc
    CDIAGocc = array2_init(occAO)
    call array2_matmul(Ctmp,Uocc,CDIAGocc,'n','n',1.0E0_realk,0.0E0_realk)
    call array2_free(Ctmp)





    ! ***************************************************************************************
    ! *                                  VIRTUAL ORBITAL SPACE                              *
    ! ***************************************************************************************

    ! For the virtual space we do exactly the same with the virt-virt blocks of the Fock matrix
    ! and the virtual local MO coefficients (MyFragment%Cv). (Thus, comments are omitted below)


    Uvirt = array2_init(virtvirt)

    call mem_alloc(S,nvirt,nvirt)
    S=0.0E0_realk
    do i=1,nvirt
       S(i,i) = 1E0_realk
    end do
    call solve_eigenvalue_problem(nvirt,MyFragment%qqfock,S,EVvirt,Uvirt%val)
    call mem_dealloc(S)

    Ctmp = array2_init(virtAO,MyFragment%Cv)
    CDIAGvirt = array2_init(virtAO)
    call array2_matmul(Ctmp,Uvirt,CDIAGvirt,'n','n',1.0E0_realk,0.0E0_realk)
    call array2_free(Ctmp)



  end subroutine get_MP2_integral_transformation_matrices

  !> \brief Construct transformation matrices used in ccsd(t) integral scheme:
  !> 1. Transformation matrix between diagonal basis (AOS) and localized basis (AOS).
  !> 2. Transformation matrix from atomic basis (atomic extent) to Fock-diagonal basis (AOS).
  !> Matrices are constructed separately for occupied and virtual spaces.
  !> Note: The array2 structures are also initialized here according to the incomming dimensions.
  !> \author Janus Juul Eriksen (adapted from mp2 routine by Kasper Kristensen),changed by PE
  !> \date February 2013
  subroutine get_canonical_integral_transformation_matrices(no,nv,nb,ppfock,qqfock,&
             &Co,Cv,CDIAGocc, CDIAGvirt, Uocc, Uvirt,EVocc, EVvirt)

    implicit none

    !> nocc, nvirt, and nbasis for fragment or full molecule
    integer, intent(in) :: no, nv, nb
    !> ppfock and qqfock for fragment or full molecule
    real(realk), intent(in) :: ppfock(no,no), qqfock(nv,nv)
    !> mo coefficents for occ and virt space for fragment or full molecule
    real(realk), intent(in) :: Co(nb,no), Cv(nb,nv)
    !> Transforming from AO to occupied orbitals in the basis where the fock matrix is diagonal
    !> - indices: (aobasis,psuedo-canonical)
    real(realk),intent(inout) :: CDIAGocc(nb,no)
    !> Transforming from AO to virtual orbitals in the basis where the fock matrix is diagonal
    !> - indices: (aobasis,pseudo-canonical)
    real(realk),intent(inout) :: CDIAGvirt(nb,nv)
    !> Occupied space: Transforming between local basis and basis where the fock matrix is diagonal
    real(realk),intent(inout) :: Uocc(no,no)
    !> Virtual space: Transforming between local basis and basis where the fock matrix is diagonal
    real(realk),intent(inout) :: Uvirt(nv,nv)
    !> Eigenvalues from diagonalization of occ-occ Fock matrix block
    real(realk), dimension(no),intent(inout) :: EVocc
    !> Eigenvalues from diagonalization of virt-virt Fock matrix block
    real(realk), dimension(nv),intent(inout) :: EVvirt
    real(realk), pointer :: S(:,:)
    integer :: i!,nocc,nvirt,nbasis,i
    integer, dimension(2) :: occocc,virtvirt,occAO,virtAO
    type(array2) :: Ctmp
    

    ! Init stuff
    ! **********
    occocc   = [no,no]
    virtvirt = [nv,nv]
    occAO    = [nb,no]
    virtAO   = [nb,nv]

    ! ****************************************************************************************
    ! *                                  OCCUPIED ORBITAL SPACE                              *
    ! ****************************************************************************************

    ! 1. Solve Fock eigenvalue problem
    ! ********************************


    ! The overlap matrix is simply the unit matrix because the local basis functions are orthogonal.
    call mem_alloc(S,no,no)
    S=0.0E0_realk
    do i=1,no
       S(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem

    call solve_eigenvalue_problem(no,ppfock,S,EVocc,Uocc)

    call mem_dealloc(S)

    ! Now Uocc contains the transformation coefficients between the local and diagonal basis.
    ! In particular, using capital/small letters for diagonal/local basis we have:
    ! phi_{K} = sum_{k} phi_{k} U_{k K}
    ! phi_{k} = sum_{K} phi_{K} U_{k K}


    ! 2. From AO to diagonal basis
    ! ****************************

    ! The transformation from the AO basis to the diagonal basis CDIAGocc is then:
    !
    ! CDIAGocc = CLOCALocc Uocc     (*)
    !
    ! where CLOCALocc are the basis MO coefficients in MyFragment%Co

    ! Determine CDIAGocc
    call dgemm('n','n',nb,no,no,1.0E0_realk,Co,nb,Uocc,no,0.0E0_realk,CDIAGocc,nb)

    ! ***************************************************************************************
    ! *                                  VIRTUAL ORBITAL SPACE                              *
    ! ***************************************************************************************

    ! For the virtual space we do exactly the same with the virt-virt blocks of the Fock matrix
    ! and the virtual local MO coefficients (MyFragment%Cv). (Thus, comments are omitted below)

    call mem_alloc(S,nv,nv)
    S=0.0E0_realk
    do i=1,nv
       S(i,i) = 1E0_realk
    end do
    call solve_eigenvalue_problem(nv,qqfock,S,EVvirt,Uvirt)
    call mem_dealloc(S)

    call dgemm('n','n',nb,nv,nv,1.0E0_realk,Cv,nb,Uvirt,nv,0.0E0_realk,CDIAGvirt,nb)


  end subroutine get_canonical_integral_transformation_matrices


  !> \brief Construct transformation matrices used in MP2 integral scheme for MP2 gradient
  !> using frozen core appriximation.
  !> The same coefficient matrices as in get_MP2_integral_transformation_matrices 
  !> are produced here - but only for the valence space.
  !> In addition three coefficient matrices for both core+valence space are constructed.
  !> Without frozen core approximation or for simple energy calculation WITH frozen core approx
  !> the simpler get_MP2_integral_transformation_matrices routine should be used.
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine get_MP2_integral_transformation_matrices_fc(MyFragment,CDIAGocc,CDIAGvirt, Uocc, Uvirt, &
       & LoccTALL, CDIAGoccTALL, UoccALL, EVocc, EVvirt)


    implicit none
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    !> Transforming from AO to valence occupied orbitals in Fock-diagonal basis
    type(array2),intent(inout) :: CDIAGocc
    !> Transforming from AO to virtual orbitals in Fock-diagonal basis
    type(array2),intent(inout) :: CDIAGvirt
    !> Occupied valence space: Transforming between local basis and Fock-diagonal basis
    type(array2),intent(inout) :: Uocc
    !> Virtual space: Transforming between local basis and Fock-diagonal basis
    type(array2),intent(inout) :: Uvirt
    !> (Transposed) local occupied MO coefficients for BOTH core and valence
    type(array2),intent(inout) :: LoccTALL
    !> Same as CDIAGocc but both for core+valence 
    !> NOTE! Here and nowhere else we put valence before core to
    !> have the most efficient storage in MP2_integrals_and_amplitudes_workhorse
    type(array2),intent(inout) :: CDIAGoccTALL
    !> Same as Uocc but for both core AND valence orbitals
    type(array2),intent(inout) :: UoccALL
    !> Eigenvalues from diagonalization of valence-valence Fock matrix block 
    real(realk), dimension(MyFragment%noccAOS),intent(inout) :: EVocc
    !> Eigenvalues from diagonalization of virt-virt Fock matrix block
    real(realk), dimension(MyFragment%nvirtAOS),intent(inout) :: EVvirt
    real(realk), pointer :: S(:,:), EVcore(:), Ucore(:,:)
    integer :: nval,nvirt,nbasis,i,nocctot,j,ncore
    integer, dimension(2) :: virtvirt,occAO,virtAO,LoccT_dims
    type(array2) :: LoccALL,Ctmp

    ! Only use for frozen core
    if(.not. DECinfo%frozencore) then
       call lsquit('get_MP2_integral_transformation_matrices_fc works only for frozen core!',-1)
    end if
 

    ! Init stuff
    ! **********
    ncore = MyFragment%ncore   ! number of core orbitals
    nval = MyFragment%noccAOS      ! (occ AOS valence size)
    nvirt = MyFragment%nvirtAOS   ! (virt AOS size)
    nbasis = MyFragment%nbasis        ! (atomic extent)
    nocctot = MyFragment%nocctot   ! occ AOS core+valence
    virtvirt = [nvirt,nvirt]
    occAO = [nbasis,nocctot]
    virtAO = [nbasis,nvirt]


    if(ncore==0) then
       call lsquit('get_MP2_integral_transformation_matrices_fc: Zero core orbitals!',-1)
    end if


    ! ****************************************************************************************
    ! *                                  OCCUPIED ORBITAL SPACE                              *
    ! ****************************************************************************************

    ! Solve Fock eigenvalue problem for core-core block
    ! *************************************************

    ! The overlap matrix is simply the unit matrix because the local basis functions are orthogonal.
    ! Note: This is not strictly true because we used fitted MO coefficients, however
    !       this approximation which is (hopefully) much less severe than the 
    !       approximation associated with the FOT.
    call mem_alloc(S,ncore,ncore)
    S=0.0E0_realk
    do i=1,ncore
       S(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem
    call mem_alloc(EVcore,ncore)
    call mem_alloc(Ucore,ncore,ncore)
    call solve_eigenvalue_problem(ncore,MyFragment%ccfock,S,EVcore,Ucore)
    call mem_dealloc(S)
    call mem_dealloc(EVcore)


    ! Solve Fock eigenvalue problem for valence-valence block
    ! *******************************************************

    call mem_alloc(S,nval,nval)
    S=0.0E0_realk
    do i=1,nval
       S(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem
    Uocc = array2_init([nval,nval])
    call solve_eigenvalue_problem(nval,MyFragment%ppfock,S,EVocc,Uocc%val)
    call mem_dealloc(S)


    ! Construct U matrix for BOTH core and valence
    ! --------------------------------------------

    ! U has the structure:
    ! 
    !         Ucore    0
    ! U =       0     Uval
    UoccALL = array2_init([nocctot,nocctot])

    ! Core
    do j=1,ncore
       do i=1,ncore
          UoccALL%val(i,j) = Ucore(i,j)
       end do
    end do
    call mem_dealloc(Ucore)

    ! Valence
    do j=1,nval
       do i=1,nval
          UoccALL%val(i+ncore,j+ncore) = Uocc%val(i,j)
       end do
    end do


    ! Local MO coefficient matrix (core+valence)
    LoccALL = array2_init(occAO)
    do j=1,ncore     ! Put core orbitals into LoccALL
       LoccALL%val(:,j) = MyFragment%coreMO(:,j)
    end do
    do j=1,nval      ! Put valence orbitals into LoccALL
       LoccALL%val(:,j+ncore) = MyFragment%Co(:,j)
    end do


    ! From AO to diagonal basis
    ! *************************

    ! The transformation from the AO basis to the diagonal basis is:
    !
    ! CDIAG = CLOCAL U           (*)

    ! Determine CDIAG for occupied space (core+valence) and store in Ctmp
    Ctmp = array2_init(occAO)
    call array2_matmul(LoccALL,UoccALL,Ctmp,'n','n',1.0E0_realk,0.0E0_realk)

    ! Put all CDIAG coefficients into transposed matrix
    ! --------------------------------------------------
    ! SPECIAL NOTE!!! Here, and nowhere else, we put valence before core!
    ! This is done to have the optimal ordering with respect to consecutive memory storage
    ! inside MP2_integrals_and_amplitudes_workhorse.
    CDIAGoccTALL = array2_init([nocctot,nbasis])

    ! Valence
    do j=1,nbasis
       do i=1,nval
          CDIAGoccTALL%val(i,j) = Ctmp%val(j,i+ncore)
       end do
    end do

    ! Core
    do j=1,nbasis
       do i=1,ncore
          CDIAGoccTALL%val(i+nval,j) = Ctmp%val(j,i)
       end do
    end do

    ! Also put valence coefficients into non-transposed matrix
    ! --------------------------------------------------------
    CDIAGocc = array2_init([nbasis,nval])
    do i=1,nval
       CDIAGocc%val(:,i) = Ctmp%val(:,i+ncore)
    end do
    call array2_free(Ctmp)


    ! Copy local occupied MOs for both core and valence into transposed array
    LoccT_dims = [nocctot,nbasis]
    LoccTALL = array2_init(LoccT_dims)

    ! Put core orbitals into LoccTALL
    do i=1,nbasis
       do j=1,nocctot
          LoccTALL%val(j,i) = LoccALL%val(i,j)
       end do
    end do
    call array2_free(LoccALL)




    ! ***************************************************************************************
    ! *                                  VIRTUAL ORBITAL SPACE                              *
    ! ***************************************************************************************

    ! For the virtual space we do the same with the virt-virt blocks of the Fock matrix.


    Uvirt = array2_init(virtvirt)

    call mem_alloc(S,nvirt,nvirt)
    S=0.0E0_realk
    do i=1,nvirt
       S(i,i) = 1E0_realk
    end do
    call solve_eigenvalue_problem(nvirt,MyFragment%qqfock,S,EVvirt,Uvirt%val)
    call mem_dealloc(S)

    Ctmp = array2_init(virtAO,MyFragment%Cv)
    CDIAGvirt = array2_init(virtAO)
    call array2_matmul(Ctmp,Uvirt,CDIAGvirt,'n','n',1.0E0_realk,0.0E0_realk)
    call array2_free(Ctmp)



  end subroutine get_MP2_integral_transformation_matrices_fc





end module array2_simple_operations
