

MODULE lattice_vectors
use precision
use matrix_module, only: matrix
use typedeftype, only: lssetting
use molecule_typetype, only: MOLECULEINFO
use molecule_type, only: init_Moleculeinfo, copy_atom
use matrix_operations, only: mat_free
use lattice_type

CONTAINS


SUBROUTINE pbc_setup_default(latt_config)
IMPLICIT NONE
TYPE(lvec_list_t),INTENT(INOUT) :: latt_config

latt_config%COMP_PBC=.false.
latt_config%compare_elmnts=.false.
latt_config%setup_pbclatt=.false.

END SUBROUTINE pbc_setup_default

SUBROUTINE READ_LATT_VECTORS(LUPRI,LUINFO,ll)
IMPLICIT NONE
INTEGER, INTENT(IN) :: LUPRI, LUINFO
TYPE(lvec_list_t), INTENT(INOUT) :: ll
CHARACTER(len=80) :: TEMPLINE
INTEGER :: IPOS,IPOS2
CHARACTER(len=10) :: activedim

ll%ldef%is_active=.true.
READ(LUINFO, '(a80)') TEMPLINE
IPOS=INDEX(TEMPLINE,'a1')
IF(IPOS .gt. 0) THEN
  IPOS2=INDEX(TEMPLINE(IPOS:),'=')
  IF(IPOS2 .eq. 0 .or. IPOS2 .GT. 7) THEN
    write(LUPRI,*) 'Incorrect usage of lattice vectors'
    write(LUPRI,*) 'Format is a1= a11 a12 a13'
    call LSQUIT('Incorrect input for lattice vectors',lupri)
  ENDIF
  READ(TEMPLINE(IPOS+IPOS2:80),*) ll%ldef%avec(1,1),ll%ldef%avec(2,1),ll%ldef%avec(3,1), activedim
if(activedim=='inactive') ll%ldef%is_active(1)= .false.
ELSE
    call LSQUIT('Incorrect input for lattice vectors',lupri)
ENDIF
READ(LUINFO, '(a80)') TEMPLINE
IPOS=INDEX(TEMPLINE,'a2')
IF(IPOS .gt. 0) THEN
  IPOS2=INDEX(TEMPLINE(IPOS:),'=')
  IF(IPOS2 .eq. 0 .or. IPOS2 .GT. 7) THEN
    write(LUPRI,*) 'Incorrect usage of lattice vectors'
    write(LUPRI,*) 'Format is a2= a21 a22 a23'
    call LSQUIT('Incorrect input for lattice vectors',lupri)
  ENDIF
  READ(TEMPLINE(IPOS+IPOS2:80),*) ll%ldef%avec(1,2),ll%ldef%avec(2,2),ll%ldef%avec(3,2), activedim
 if(activedim=='inactive') ll%ldef%is_active(2)=.false.
ENDIF
READ(LUINFO, '(a80)') TEMPLINE
IPOS=INDEX(TEMPLINE,'a3')
IF(IPOS .gt. 0) THEN
  IPOS2=INDEX(TEMPLINE(IPOS:),'=')
  IF(IPOS2 .eq. 0 .or. IPOS2 .GT. 7) THEN
    write(LUPRI,*) 'Incorrect usage of lattice vectors'
    write(LUPRI,*) 'Format is a3= a31 a32 a33'
    call LSQUIT('Incorrect input for lattice vectors',lupri)
  ENDIF
 READ(TEMPLINE(IPOS+IPOS2:80),*) ll%ldef%avec(1,3),ll%ldef%avec(2,3),ll%ldef%avec(3,3), activedim
 if(activedim=='inactive') ll%ldef%is_active(3)= .false.
ENDIF
!call write_matrix(ll%ldef%avec,3,3)

END SUBROUTINE READ_LATT_VECTORS

Logical function ifpbc_active(n,lupri)
INTEGER, INTENT(IN) :: n,lupri

  if(n .lt. 1 .and. n .gt. 3) call LSQUIT('if_active, n must be: 1,2 or 3',lupri)
  ifpbc_active=pbc_control%ldef%is_active(n)

END function ifpbc_active

SUBROUTINE translate_atom(setting,lattice_cell,iao,natoms)
  IMPLICIT NONE
  !LOGICAL, INTENT(IN) :: atom2trans(4)
  TYPE(LSSETTING), INTENT(INOUT) ::setting
  INTEGER, INTENT(IN) ::  natoms
  INTEGER :: iatom
  TYPE(lattice_cell_info_t), INTENT(IN)  :: lattice_cell
  INTEGER,INTENT(IN) :: iao
  REAL(realk) :: x,y, z

      DO iatom=1,natoms
      
         x=lattice_cell%atom(iatom)%center(1)
         y=lattice_cell%atom(iatom)%center(2)
         z=lattice_cell%atom(iatom)%center(3)
         setting%molecule(iao)%p%atom(iatom)%center(1)=x !&
         setting%molecule(iao)%p%atom(iatom)%center(2)=y !&
         setting%molecule(iao)%p%atom(iatom)%center(3)=z !&
         !setting%molecule(iao)%p%atom(iatom)%center=&
          !& lattice_cell%atom(iatom)%center

      ENDDO

END SUBROUTINE translate_atom

SUBROUTINE set_lattice_cells(lattice_cell,num_latvectors,molecule,ll,lupri)
  implicit none
  ! local variables
  INTEGER :: IATOM, index
  INTEGER, INTENT(IN) :: num_latvectors,lupri
  TYPE(MOLECULEINFO), INTENT(IN) :: molecule
  TYPE(moleculeinfo), INTENT(INOUT), DIMENSION(num_latvectors) :: lattice_cell
  TYPE(lvec_list_t),intent(in)  :: ll
  CHARACTER(len=22) :: mollabel
   
!  Allocate(lattice_cell(num_latvectors))
  write(lupri,*) 'Number of atoms ', molecule%natoms
  !call build_lvec_list(ll)
!     write(lupri,*) 'before loop'
  DO index=1,num_latvectors
!  DO index=-ll%max_layer,max_layer
!  Allocate(lattice_cell(index)%atom(molecule%natoms))
     write(mollabel,'(A12,I9)') 'PBC-Molecule',index
     call init_Moleculeinfo(lattice_cell(index),molecule%natoms,mollabel)


!     write(*,*) 'before copy_atom'
     DO iatom=1,molecule%natoms
      call copy_atom(molecule,iatom,lattice_cell(index),iatom,6)
     ENDDO
!     write(*,*) 'after copy_atom',molecule%natoms
     DO IATOM=1,molecule%natoms
!     write(*,*) 'inside loop copy_atom'
        lattice_cell(index)%atom(IATOM)%CENTER(1)= &
        &  molecule%atom(IATOM)%CENTER(1)-ll%lvec(index)%std_coord(1)
!     write(*,*) 'inside loop copy_atom'

        lattice_cell(index)%atom(IATOM)%CENTER(2)= &
        &  molecule%atom(IATOM)%CENTER(2)-ll%lvec(index)%std_coord(2)

        lattice_cell(index)%atom(IATOM)%CENTER(3)= &
        &  molecule%atom(IATOM)%CENTER(3)-ll%lvec(index)%std_coord(3)
!        write(*,*) 'x', lattice_cell(index)%atom(IATOM)%CENTER(1)
!        write(*,*) 'y', lattice_cell(index)%atom(IATOM)%CENTER(2)
!        write(*,*) 'z', lattice_cell(index)%atom(IATOM)%CENTER(3)
     ENDDO
!     write(*,*) 'after loop copy_atom'
  ENDDO
!     write(*,*) 'after second loop copy_atom'
  write(lupri,*) 'Finished set_lattice_cells'
!  stop

END SUBROUTINE set_lattice_cells

SUBROUTINE newpbc_get_molecules(setting,molecule)
  implicit none
!  #include <priunit.h>
!  #include <inforb.h>
!  #include <pbc.h>
  ! input and output arguments

  ! local variables
  TYPE(LSSETTING), INTENT(INOUT)   :: SETTING
  INTEGER :: IAO,IATOM
  TYPE(MOLECULEINFO), INTENT(IN) :: molecule


    DO iao=1,4
     NULLIFY(SETTING%MOLECULE(IAO)%p)
!CAREFUL Johannes husk å deallokere denne (lag free_pbc_molecules)
     ALLOCATE(SETTING%MOLECULE(IAO)%p)
!  write(lupri,*) 'Test overlap matrix'
     call init_Moleculeinfo(setting%molecule(iao)%p,molecule%natoms,&
    &                       'PBC-Molecule          ')
!  write(lupri,*) 'Test overlap matrix'
     DO IATOM=1,molecule%natoms
        call copy_atom(molecule,iatom,setting%molecule(iao)%p,iatom,6)
     ENDDO
!     NULLIFY(SETTING%MOLECULE(IAO)%p)
     setting%molbuild(iao)=.true.
     setting%fragment(iao)%p => setting%molecule(iao)%p

   ENDDO

END SUBROUTINE newpbc_get_molecules

SUBROUTINE set_refcell(refcell,molecule)
  implicit none
  INTEGER :: IATOM 
  TYPE(MOLECULEINFO), INTENT(IN) :: molecule
  TYPE(MOLECULEINFO), INTENT(INOUT) :: refcell

  call init_Moleculeinfo(refcell,molecule%natoms,&
    &                       'PBC-Reference-cell    ')

     DO IATOM=1,molecule%natoms
        call copy_atom(molecule,iatom,refcell,iatom,6)
     ENDDO


END SUBROUTINE set_refcell


!This subroutine will setup the lattice vectors
SUBROUTINE build_lvec_list(ll,nbast)
  implicit none
  INTEGER, intent(in) ::nbast
  type(lvec_list_t), intent(inout) :: ll
  INTEGER:: l1, l2,l3, fdim(3),alstat,index

  ! invent some basic lattice data
!  ll%ldef%avec(1:3,1) = (/ 9.2822, 0.0, 0.0 /)
!  ll%ldef%avec(1:3,2) = (/ 0.0, 10.0, 0.0 /)
!  ll%ldef%avec(1:3,3) = (/ 0.0, 0.0, 10.0 /)
!  ll%ldef%is_active(1:3) = (/ .true., .false., .false. /)

  ! generate a very crude list of lattice vectors corresponding to a
  ! block of cells
!  write(*,*) 'avec'
!  call write_matrix(ll%ldef%avec,3,3)
  fdim(1:3) = 0
  if (ll%ldef%is_active(1)) fdim(1) = 1
  if (ll%ldef%is_active(2)) fdim(2) = 1
  if (ll%ldef%is_active(3)) fdim(3) = 1


!  max_layer = 1
  ll%num_entries = (2*ll%max_layer*fdim(1)+1)*(2*ll%max_layer*fdim(2)+1)*(2*ll%max_layer*fdim(3)+1) 
  allocate(ll%lvec(ll%num_entries), STAT=alstat)
!  allocate(ll%lvec(-ll%max_layer:ll%max_layer), STAT=alstat)
!  if (alstat .eq. ???) then
!     quit('mem alloc failed')
!  end if

  index = 1
!  index = -ll%max_layer

  do l3 = -ll%max_layer*fdim(3), ll%max_layer*fdim(3)
  do l2 = -ll%max_layer*fdim(2), ll%max_layer*fdim(2)
  do l1 = -ll%max_layer*fdim(1), ll%max_layer*fdim(1)
     ll%lvec(index)%lat_coord(1:3) = (/ real(l1), real(l2), real(l3) /)
     call latt_2_std_coord(ll%lvec(index)%lat_coord,ll%lvec(index)%std_coord,ll%ldef%avec)  
     allocate(ll%lvec(index)%fck_vec(nbast*nbast))
     allocate(ll%lvec(index)%fck_mat(nbast,nbast))
     allocate(ll%lvec(index)%Sl_vec(nbast*nbast))
     allocate(ll%lvec(index)%Sl_mat(nbast,nbast))
     allocate(ll%lvec(index)%d_vec(nbast*nbast))
     allocate(ll%lvec(index)%d_mat(nbast,nbast))
     ll%lvec(index)%fck_vec=0 
     ll%lvec(index)%fck_mat=0 
     ll%lvec(index)%Sl_vec=0 
     ll%lvec(index)%Sl_mat=0 
     ll%lvec(index)%is_redundant=.false.
     index = index + 1
  end do
  end do
  end do
  ! remember to deallocate ll%lvec somewhere in the code!

  ll%opdat(1)%basename='ovl'
  ll%opdat(2)%basename='kin'
  ll%opdat(3)%basename='nucn'
  ll%opdat(4)%basename='coul'
  ll%opdat(5)%basename='exch'
  ll%opdat(6)%basename='xcor'
  ll%opdat(7)%basename='fock'
  ll%opdat(8)%basename='dmat'
  ll%opdat(9)%basename='nucf'


end subroutine build_lvec_list

SUBROUTINE build_nflvec_list(ll,nbast)
  implicit none
  INTEGER, intent(in) ::nbast
  type(lvec_list_t), intent(inout) :: ll
  INTEGER:: l1, l2,l3, fdim(3),alstat,index

  fdim(1:3) = 0
  if (ll%ldef%is_active(1)) fdim(1) = 1
  if (ll%ldef%is_active(2)) fdim(2) = 1
  if (ll%ldef%is_active(3)) fdim(3) = 1


!  max_layer = 1
  ll%nf_entries = (2*ll%nneighbour*fdim(1)+1)*(2*ll%nneighbour*fdim(2)+1)*(2*ll%nneighbour*fdim(3)+1) 
  allocate(ll%nflvec(ll%nf_entries), STAT=alstat)

  index=1
  do l3 = -ll%nneighbour*fdim(3), ll%nneighbour*fdim(3)
  do l2 = -ll%nneighbour*fdim(2), ll%nneighbour*fdim(2)
  do l1 = -ll%nneighbour*fdim(1), ll%nneighbour*fdim(1)
     ll%nflvec(index)%lat_coord(1:3) = (/ real(l1), real(l2), real(l3) /)
     call latt_2_std_coord(ll%nflvec(index)%lat_coord,ll%nflvec(index)%std_coord,ll%ldef%avec)  
     allocate(ll%nflvec(index)%fck_vec(nbast*nbast))
     allocate(ll%nflvec(index)%fck_mat(nbast,nbast))
     allocate(ll%nflvec(index)%d_vec(nbast*nbast))
     allocate(ll%nflvec(index)%d_mat(nbast,nbast))
     ll%nflvec(index)%fck_vec=0 
     ll%nflvec(index)%fck_mat=0 
     index = index + 1
  end do
  end do
  end do



END SUBROUTINE build_nflvec_list

subroutine reset_lvec_list(ll)
  implicit none

  type(lvec_list_t), intent(inout) :: ll

  ll%ldef%is_active(1:3) = .false.
  ll%ldef%avec(1:3,1:3) = 0.0
  ll%num_entries = 0
  deallocate(ll%lvec)
  nullify(ll%lvec)

END SUBROUTINE reset_lvec_list

SUBROUTINE free_lvec_list(lvec_list)
implicit none
TYPE(lvec_list_t),intent(inout) :: lvec_list
integer :: i,j
integer :: x,y,z

!call free_lvec_data(lvec_list%lvec,size(lvec_list%lvec))
!call free_lvec_data(lvec_list%nflvec,size(lvec_list%nflvec))
 !do i=1,size(lvec_list%lvec)
 !   call find_latt_vectors(i,x,y,x,fdim,lvec_list)
 deallocate(lvec_list%lvec)
    

END SUBROUTINE free_lvec_list

SUBROUTINE free_lvec_data(lvec_data,numvecs)
implicit none
INTEGER,INTENT(IN) :: numvecs
TYPE(lvec_data_t),intent(inout) :: lvec_data(numvecs)
!local variables
integer :: i,j

do i=1,numvecs
 do j=1,MaxPBCOpTypes
!  if(associated(lvec_data(i)%oper(j))) then
	call mat_free(lvec_data(i)%oper(j))
!  endif
 enddo
enddo

END SUBROUTINE free_lvec_data

!Transforms from lattice coordinates to standard coordinates
SUBROUTINE latt_2_std_coord(latcoord,stdcoord,latvec)
  implicit none
  ! input and output arguments
  real(realk), intent(in)    :: latvec(3,3),latcoord(3)
  real(realk), intent(inout)   :: stdcoord(3)
  ! local variables
  integer             :: ci, vi

!  if (.not. lat_data%latvec_ok) then
 !    call quit('Call to PBC_LAT2STD_COORD, but no lattice vectors have been supplied to the PBC_DATA module.')
 ! end if
!     write(*,*) 'll%avec', latvec
!     write(*,*) 'latcoord', latcoord

  do ci = 1,3
     stdcoord(ci) = 0.0
!     stdcoord(ci) = latvec(ci,1)+ latvec(ci,2)+latvec(ci,3)
     do vi = 1,3
        stdcoord(ci) = stdcoord(ci) + latcoord(vi) * latvec(ci,vi)
!        write(*,*) 'latvec(ci,vi)', latvec(vi,ci)
     end do
   !write(*,*) 'stdcoord', stdcoord(ci)
  end do
   !write(*,*) 'stdcoord', stdcoord
END SUBROUTINE latt_2_std_coord


!Find the lattice vectors from the dummy lattice index ll.
SUBROUTINE find_latt_vectors(ll,l1,l2,l3,fdim,latt)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: fdim(3)
INTEGER, INTENT(IN) :: ll
INTEGER, INTENT(INOUT) :: l1,l2,l3
TYPE(lvec_list_t), INTENT(IN) :: latt
!!!!!
!INTEGER :: vec1,vec2,vec3, idex

!fdim(1:3) = 0
!if (latt%ldef%is_active(1)) fdim(1) = 1
!if (latt%ldef%is_active(2)) fdim(2) = 1
!if (latt%ldef%is_active(3)) fdim(3) = 1
l1=int(latt%lvec(ll)%lat_coord(1))
l2=int(latt%lvec(ll)%lat_coord(2))
l3=int(latt%lvec(ll)%lat_coord(3))


!idex=0
!DO vec3=-fdim(3)*latt%max_layer,fdim(3)*latt%max_layer
! DO vec2=-fdim(2)*latt%max_layer,fdim(2)*latt%max_layer
!  DO vec1=-fdim(1)*latt%max_layer,fdim(1)*latt%max_layer
!    idex=idex+1
!    IF(idex .eq. ll ) THEN
!      l1=vec1
!      l2=vec2
!      l3=vec3
!      EXIT
!    ENDIF
!  ENDDO
! ENDDO
!ENDDO

END SUBROUTINE find_latt_vectors

SUBROUTINE pbcstruct_get_active_dims(fdim)
implicit none
integer, intent(inout) :: fdim(3)

fdim(:) = 0
  if (pbc_control%ldef%is_active(1)) fdim(1) = 1
  if (pbc_control%ldef%is_active(2)) fdim(2) = 1
  if (pbc_control%ldef%is_active(3)) fdim(3) = 1
  
END SUBROUTINE pbcstruct_get_active_dims

!Transforms from the lattice vectors l1, l2 and l3 to the lattice index ll
SUBROUTINE find_latt_index(ll,l1,l2,l3,fdim,latt,maxlayer)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: fdim(3)
INTEGER, INTENT(IN) :: l1,l2,l3
INTEGER, INTENT(IN) :: maxlayer
INTEGER, INTENT(OUT) :: ll
TYPE(lvec_list_t),INTENT(IN) :: latt
!!!!!
INTEGER :: nx,ny,nz

fdim(1:3) = 0
  if (latt%ldef%is_active(1)) fdim(1) = 1
  if (latt%ldef%is_active(2)) fdim(2) = 1
  if (latt%ldef%is_active(3)) fdim(3) = 1

  nx=maxlayer*fdim(1)
  ny=maxlayer*fdim(2)
  nz=maxlayer*fdim(3)

  ll=l1+nx+1+2*(ny+l2)*nx+ny+l2+2*(nz+l3)*(nx+ny)+4*(nz+l3)*nx*ny+nz+l3


!idex=0
!DO vec1=-fdim(1)*latt%max_layer,fdim(1)*latt%max_layer
! DO vec2=-fdim(2)*latt%max_layer,fdim(2)*latt%max_layer
!  DO vec3=-fdim(3)*latt%max_layer,fdim(3)*latt%max_layer
!    idex=idex+1
!    IF(vec1 /= l1) CYCLE
!    IF(vec2 /= l2) CYCLE
!    IF(vec3 /= l3) CYCLE
!    ll=idex
!  ENDDO
! ENDDO
!ENDDO





END SUBROUTINE find_latt_index


END MODULE lattice_vectors

