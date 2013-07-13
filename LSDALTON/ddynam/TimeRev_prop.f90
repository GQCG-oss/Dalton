!==============================!
! Time-reversible propagation  !
!==============================!
Module TimeRev_propagation
! Time reversible propagation of electronic degrees
! of freedom
! A.M.N.Niklasson et al. PRL 97, 123001(2006)
Use precision
Use ls_util
Use matrix_module
Use matrix_operations
Contains
!==============!
! Propagation  !
!==============!
! Propagates the density matrix
Subroutine Propagation(N,nbas,StepNum,D,Darr,Daux,Start_propagation)
!
Implicit none
Integer :: N ! Number of fitted values
Integer :: nbas,StepNum ! Order of D matrix, step number 
Type(matrix), intent(inout) :: D
Real(realk), dimension(nbas,nbas) :: Dunpacked
Real(realk) Daux(:,:,:) ! Array of auxiliary density matrices
Real(realk) Darr(:,:,:) ! Array of converged density matrices
Logical :: Start_propagation
Write(*,*)'StepNum=',StepNum
! Unpacking D
Call Mat_to_full(D,1E0_realk,Dunpacked)
Write(*,*) 'Dconv=', Dunpacked(1,1)
! Filling in Daux
If (StepNum .LT. N) then
   Daux(StepNum+1,:,:) = Dunpacked
Endif
! Filling in Darr
If (N .EQ. 6) then
   If ((StepNum .LT. N) .AND. (StepNum .GE. 1)) then
      Darr(StepNum,:,:) = Dunpacked
   Else
      Darr = eoshift(Darr,Shift=1,boundary=Dunpacked,dim=1)
   Endif
Endif
! Propagation starts
   ! Two fitted values
If (N .EQ. 2) then
   If (StepNum .GE. (N-1)) then
   Start_propagation = .TRUE.
      Dunpacked = 2*Dunpacked - Daux(N-1,:,:)
      Daux = eoshift(Daux,Shift=1,boundary=Dunpacked,dim=1)
      Write(*,*) 'Dprop=', Dunpacked(1,1)
      ! Convert Dunpacked to a matrix type
      Call mat_set_from_full(Dunpacked,1E0_realk,D)
   Endif
Endif
   ! Six fitted values
If (N .EQ. 6) then
   If (StepNum .GE. N) then
      Dunpacked = -((1.0D0/13.0D0)*(30.0D0*(Darr(N-1,:,:) + Darr(N-5,:,:))-&
      & 3.0D0*(Darr(N-2,:,:)+Darr(N-4,:,:))-28.0D0*Darr(N-3,:,:))&
      & -Daux(N-5,:,:))
      Daux = eoshift(Daux,Shift=1,boundary=Dunpacked,dim=1)
      Write(*,*) 'Dprop=', Dunpacked(1,1)
      ! Convert Dunpacked to a matrix type
      Call mat_set_from_full(Dunpacked,1E0_realk,D)
   Endif
Endif
!
End subroutine Propagation
!
End module TimeRev_propagation
