#ifdef MOD_UNRELEASED
module pbcffdata
  use precision
  use TYPEDEF
  implicit none
  public

  integer, parameter :: min_qlmax = 3, max_qlmax = 21
  integer, parameter :: min_TWlmax = 5, max_TWlmax = 50

  type ffdata_t
     logical :: reset_T_diverg
     logical :: square_intermediates
     integer :: q_lmax     ! Lmax for multipole moments
     integer :: TW_lmax    ! Lmax for T and W tensors
     real(realk), pointer :: cell_locmom(:)
     real(realk), pointer :: Tlattice(:,:)
  end type ffdata_t

  type(ffdata_t), save :: ffdata

  contains 

  SUBROUTINE READ_pbc_tlattice(tlattice,lmax,mattxt,lupri)
  IMPLICIT NONE
  INTEGER,intent(IN) :: lmax,lupri
  real(realk),intent(Inout) ::tlattice((lmax+1)**2,(lmax+1)**2)
  character(len=18),intent(in) :: mattxt
  !LOCAL
  INTEGER :: i,j,dimin,dimtest
  INTEGER :: iunit

  dimin=(lmax+1)**2
  iunit=-1

  write(*,*) mattxt
  !write(mattxt,'(A18)') 'Tlatticetensor.dat'
  CALL lsOPEN(IUNIT,mattxt,'OLD','UNFORMATTED')
  read(iunit) dimtest!,nil
  if(dimtest .ne. dimin) then
    write(*,*) 'ERROR: dims do not match: input dim=',dimin,'read dim=',dimtest
    call lsquit('In read_pbc_tlattice: dims do not match',lupri)
  endif
 ! DO j=1,dimin
 !   read(iunit) (tlattice(j,i),i=1,dimin)
 ! ENDDO

  DO j=1,dimin
     read(iunit) (tlattice(j,i),i=1,dimin)
  enddo

  CALL lsCLOSE(IUNIT,'KEEP')
  
  END SUBROUTINE READ_pbc_tlattice

!  SUBROUTINE pbc_mmcontrib(tlat,lmax)
!  use precision
!  use TYPEDEF
!  implicit none
!  INTEGER,intent(IN) :: lmax
!  real(realk),intent(Inout) ::tlat((lmax+1)**2,(lmax+1)**2)
!
!
!
!
!  END SUBROUTINE pbc_mmcontrib
end module pbcffdata

#endif
