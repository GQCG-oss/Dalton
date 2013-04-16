!> @file
!> Simple operations for RI arrays (not used now / or any more)
!> \author Marcin Ziolkowski
module ri_simple_operations

  use precision
  use memory_handling
  use typedeftype
  use files!,only:lsopen,lsclose
  use matrix_module!, only:matrix
  use matrix_operations
  use BUILDAOBATCH!, only: shell_number,orbital_degeneracy
  use dec_typedef_module
  use integralinterfaceMOD


  !> Number of an array
  integer :: RINumber=0
  !> Number of created arrays
  integer :: CreatedRI=0
  !> Number of destroyed arrays
  integer :: DestroyedRI=0
  !> File unit counter
  integer :: RI_FileUnit=900


contains

  !> \brief Construct RI array
  function ri_init(dims) result(res)

    implicit none
    type(ri) :: res
    integer, dimension(3), intent(in) :: dims
    logical :: file_exist
    integer :: i

    do i=1,3
       res%dims(i)=dims(i)
    end do

    if(associated(res%val)) res%val => null()
    call mem_alloc(res%val,dims(1),dims(2),dims(3) )
    res%val = 0.0E0_realk

    RINumber = RINumber + 1
    CreatedRI = CreatedRI + 1
    RI_fileunit = RI_fileunit + 1

    write(res%FileName,'("ri_array_",i9.9,".data")') RINumber
    res%FUnit = RI_fileunit

    inquire(file=res%FileName,exist=file_exist)
    if(file_exist) then
       call lsopen(res%FUnit,res%FileName,'OLD','UNFORMATTED')
       call lsclose(res%FUnit,'DELETE')
    end if

    call lsopen(res%FUnit,res%FileName,'NEW','UNFORMATTED')
    call lsclose(res%FUnit,'KEEP')

    return
  end function ri_init

  !> \brief Destroy RI type
  subroutine ri_free(array)

    implicit none
    type(ri), intent(inout) :: array

    if(associated(array%val)) then
       call mem_dealloc(array%val)
       array%val => null()
    end if

    DestroyedRI = DestroyedRI + 1

    call lsopen(array%FUnit,array%FileName,'OLD','UNFORMATTED')
    call lsclose(array%FUnit,'DELETE')

    return
  end subroutine ri_free

  !> \brief RI type statistics
  subroutine ri_print_statistics(output)

    implicit none
    integer, intent(in) :: output

    write(DECinfo%output,'(/,a)')    ' ------------------------- '
    write(DECinfo%output,'(a)')      '       RI statistics       '
    write(DECinfo%output,'(a)')      ' ------------------------- '
    write(DECinfo%output,'(a,i4)')   'Number of created RI objs.   : ',CreatedRI
    write(DECinfo%output,'(a,i4)')   'Number of destroyed RI objs. : ',DestroyedRI
    write(DECinfo%output,'(a,i4)')   'Max. used unit nuber         : ',RI_fileunit
    write(DECinfo%output,'(a,i4,/)') 'Orphaned objs.               : ',CreatedRI-DestroyedRI

    return
  end subroutine ri_print_statistics

  !> \brief Print RI data
  subroutine ri_print(this,output)
    implicit none
    type(ri), intent(in) :: this
    integer, intent(in) :: output
    integer :: a,b,c


    write(DECinfo%output,'(/,a,/)') '-- data --'
    do c=1,this%dims(3)

       do a=1,this%dims(1)
          do b=1,this%dims(2)
             write(DECinfo%output,'(3i4,f16.10)') c,a,b,this%val(a,b,c)
          end do
       end do

    end do

    return
  end subroutine ri_print

  !> \brief Reset number of array parameters
  subroutine ri_reset()

    implicit none
    CreatedRI=0
    DestroyedRI=0
    RI_fileunit=900
    return
  end subroutine ri_reset

  !> \brief RI info
  subroutine ri_info(array,output)

    implicit none
    type(ri), intent(in) :: array
    real(realk) :: memory
    integer, intent(in) :: output

    write(DECinfo%output,'(/,a)')   '       RI Info      '
    write(DECinfo%output,'(a)')     '--------------------'
    write(DECinfo%output,'(a,i4)')  '    Dim1 : ',array%dims(1)
    write(DECinfo%output,'(a,i4)')  '    Dim2 : ',array%dims(2)
    write(DECinfo%output,'(a,i4)')  '    Dim3 : ',array%dims(3)
    write(DECinfo%output,'(a,i4)')  'FileUnit : ',array%FUnit
    write(DECinfo%output,'(a,a)')   'FileName : ',array%FileName

    memory=( dble(array%dims(1))*dble(array%dims(2))*dble(array%dims(3))*8 )/dble(2**20)
    write(DECinfo%output,'(a,f10.2,a)') '  Memory : ',memory,' MB'

    return
  end subroutine ri_info

  !> \brief Construct ri intermediates in AO basis
  function get_ao_ri_intermediate(mylsitem) result(l_ao)

    implicit none
    type(matrix), pointer :: rimat(:)
    type(matrix) :: v,v12_mat

    type(lsitem), intent(inout) :: mylsitem
    !    type(lssetting), intent(in) :: setting
    !    type(daltoninput), intent(in) :: dalton
    type(ri) :: ri_ao
    type(ri) :: l_ao
    integer, pointer :: AuxDegeneracy(:)
    integer :: nbas,naux
    integer :: NumAuxShell,bast
    integer :: int_output=10
    integer :: i,j,k,l,a,b
    real(realk), pointer :: v12(:,:)

    nbas = mylsitem%input%molecule%nbastREG
    naux = mylsitem%input%molecule%nbastAUX
    ri_ao = ri_init([nbas,nbas,naux])
    l_ao = ri_init([nbas,nbas,naux])

    ! get the number of auxiliary shells and their degeneracy
    call shell_number(int_output,0,mylsitem%input, &
         mylsitem%input%basis%auxiliary,NumAuxShell)
    call mem_alloc(AuxDegeneracy,NumAuxShell)
    call orbital_degeneracy(int_output,0,mylsitem%input,&
         mylsitem%input%basis%auxiliary,NumAuxShell,AuxDegeneracy)

    bast=1
    do i=1,NumAuxShell
       call mem_alloc(rimat,AuxDegeneracy(i))
       do j=1,AuxDegeneracy(i)
          call mat_init(rimat(j),nbas,nbas)
          call mat_zero(rimat(j))
       end do
       !call ii_get_3center_eri_matbatch(output,0,mylsitem%setting,bast, &
       !                AuxDegeneracy(i),rimat)
       do j=1,AuxDegeneracy(i)
          call mat_to_full(rimat(j),1.0E0_realk,ri_ao%val(:,:,bast+j-1))
       end do

       call mem_dealloc(rimat)
       bast = bast + AuxDegeneracy(i)
    end do

    call mem_dealloc(AuxDegeneracy)

    ! two-index eri
    call mat_init(v,naux,naux)
    call mat_zero(v)
    call mat_init(v12_mat,naux,naux)
    call mat_zero(v12_mat)
    !call ii_get_2center_aux_eri(output,output,mylsitem%setting,v)
    call LowdinMatrix(v,v12_mat)
    call mat_free(v)

    call mem_alloc(v12,naux,naux)
    call mat_to_full(v12_mat,1.0E0_realk,v12)

    ! construct ri intermediates
    l_ao%val = 0.0E0_realk

    do l=1,naux
       do k=1,naux
          do a=1,nbas
             do b=1,nbas
                l_ao%val(a,b,l) = l_ao%val(a,b,l) + ri_ao%val(a,b,k) * v12(k,l)
             end do
          end do
       end do
    end do

    call ri_free(ri_ao)
    call mem_dealloc(v12)

    return
  end function get_ao_ri_intermediate


  !> \brief Lowdin decomposition
  subroutine get_lowdin(msqr,n,inpm)

    implicit none
    integer, intent(in) :: n
    integer :: ifail, lwork, i !, liwork
    !    integer, dimension(128*n) :: iwork
    real(realk), intent(in), dimension(n,n) :: inpm
    real(realk), intent(out), dimension(n,n) :: msqr
    real(realk), dimension(n,n) :: x, d
    real(realk), dimension(n) :: eval
    real(realk), dimension(128*n+64*n*n) :: work
    external dsyev
    ifail=0

    !    liwork=128*n
    lwork = 128*n+64*n*n
    x = inpm
    ifail = 0
    work = 0.0E0_realk
    eval = 0.0E0_realk

    call dsyev('v','u',n,x,n,eval,work,lwork,ifail)
    !    call dsyevd('v','u',n,x,n,eval,work,lwork,iwork,liwork,ifail)
    if(ifail /= 0) then
       stop 'error :: error in dsyevd'
    end if
    d = 0.0E0_realk

    do i = 1, n
       d(i,i) = 1.0E0_realk/sqrt(eval(i))
    end do
    msqr = matmul(matmul(x,d),transpose(x))
    return
  end subroutine get_lowdin

  !> \brief Lowdin decomposition of matrix type
  subroutine LowdinMatrix(input,output)

    implicit none

    type(matrix), intent(in) :: input
    type(matrix), intent(inout) :: output
    real(realk), pointer :: input_full(:,:), output_full(:,:)

    integer :: nrow,ncol

    nrow=input%nrow
    ncol=input%ncol

    call mem_alloc(output_full,nrow,nrow)
    call mem_alloc(input_full,nrow,nrow)
    output_full=0.0E0_realk
    input_full=0.0E0_realk
    call mat_to_full(input,1.0E0_realk,input_full)

    call get_lowdin(output_full,nrow,input_full)

    call mat_zero(output)
    call mat_set_from_full(output_full,1.0E0_realk,output)
    call mem_dealloc(output_full)
    call mem_dealloc(input_full)

    return
  end subroutine LowdinMatrix

end module ri_simple_operations


