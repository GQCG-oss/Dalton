!================================!
! Transforms internal step to    !
! Cartesian coordinates          !
!================================!
Module q_to_x_mod
!
Use memory_handling
use precision
use ls_util 
use optimization_input 
use Fundamental
use molecule_type
  use lstiming, only: lstimer
Contains
!===================!
! Back_transform    !
!===================!
Subroutine Back_transform(CSTEP,N_Cart,N_Int,lupri,optinfo)
!
Implicit none
Type(opt_setting) :: optinfo
Real(realk) :: step_len
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
Real(realk), parameter :: Conv_thresh = 1.0D-12
Real(realk) CSTEP(MXCOOR)  ! Cartesian step
Real(realk), pointer :: Ini_coord(:,:), Cart_step(:),q0(:),q(:),Int_step(:),Total_ds(:)
Real(realk), pointer :: Bs_inv(:,:),Vectors(:,:),Bs(:,:),B_mat(:,:),Del_step(:), Left_ds(:)
Real(realk) :: Int_step_norm
Logical :: Converged, Finished
Integer :: N_Cart,N_Int,lupri,i,Counter,Cycles,nCent
! Initialize
step_len = D1
CSTEP = D0
Converged = .FALSE.
Finished = .FALSE.
Counter = 1
Cycles = 0
nCent = N_Cart/3
call mem_alloc(Ini_Coord,3,nCent)
call mem_alloc(Int_step,N_Int)
call mem_alloc(Del_step,N_Cart-6)
call mem_alloc(Total_ds,N_Cart-6)
call mem_alloc(Left_ds,N_Cart-6)
call mem_alloc(q0,N_Int)
call mem_alloc(q,N_Int)
call mem_alloc(Cart_step,N_Cart)
call mem_alloc(Bs_inv,N_Cart,N_Cart-6)
call mem_alloc(Bs,N_Cart-6,N_Cart)
Call mem_alloc(Vectors,N_Int,N_Cart-6)
Call mem_alloc(B_mat,N_Int,N_Cart)
!
Cart_step = D0
Ini_Coord = optinfo%Coordinates(:,1:nCent)
Int_step = optinfo%StpInt(1:N_Int)
! Get the first inverse and vectors 
Call First_inverse(optinfo,Bs_inv,vectors,N_Cart,N_Int,lupri)
!
! Transform step from redundant to delocs
!
Call DGEMV('T',N_Int,N_Cart-6,1.0E0_realk,Vectors,N_Int,&
     & Int_step,1,0.0E0_realk,Del_step,1)
! Save step in delocalized internals
Total_ds = Del_step
Left_ds = Total_ds
! Print delocalized step
If (optinfo%IPrint .GE. 12) then
  call lsheader(lupri,'Step in delocalized internals:')
  call output(Total_ds,1,1,1,N_Cart-6,1,N_Cart-6,1,LUPRI)
Endif
!
Do while (Finished .EQV. .False.)
  If (Counter .GE. 10) call lsquit('Unable to transform the step!',lupri)
  If (optinfo%New_stepping) then   ! High-order step
     Call Poly_stepping(Bs_inv,Vectors,Cart_step,Del_step,N_Cart,N_Int,N_Cart-6,Conv_thresh,&
    & Converged,step_len,lupri,optinfo) 
  Else ! Iterative back-transformation
     Call IBT(Bs_Inv,N_Cart,N_Int,N_Cart-6,optinfo,Cart_step,Del_step, &
          & Vectors,q0,lupri,Converged,Conv_thresh)
  Endif
  !
  If (.NOT. Converged) then 
     If (optinfo%IBT) call lsquit('IBT failed to transform the step!',lupri)
     Call lsheader(lupri,'BT not converged, reducing step')
     Counter = Counter + 1
     step_len = step_len*0.5E0_realk
     Del_step = Del_step*step_len
  Else  ! Converged
     Cycles = Cycles + 1
     CSTEP(1:N_Cart) = CSTEP(1:N_Cart) + Cart_step
     Do i = 1, N_Cart/3
        optinfo%Coordinates(:,i) = optinfo%Coordinates(:,i) + Cart_step(i*3-2:i*3)
     Enddo
     Converged = .FALSE.
     If (step_len .GT. 0.6) then
        Finished = .TRUE.
     Else
        step_len = D1
        Call lsheader(lupri,'Next cycle of BT')
        Left_ds = Left_ds - Del_step
        Del_step = Left_ds
        ! Get new Bs and Bs_inv using old vectors
        Call B_matrix(N_int,N_cart,2,B_mat,optinfo)
        If (optinfo%IPrint .GE. 12) then
           call lsheader(lupri,'Wilson B matrix [B(ij) = dq(i)/dx(j)] from C++')
           call output(B_mat,1,N_Int,1,N_Cart,N_Int,N_Cart,1,LUPRI)
        Endif
        Call DGEMM('T','N',N_Cart-6,N_Cart,N_Int,1.0E0_realk,&
             & Vectors,N_Int,B_mat,N_Int,0.0E0_realk,Bs,N_Cart-6)
        Call Get_B_inv(Bs,Bs_Inv,N_Cart-6,N_Cart)
        ! Some printout
        If (optinfo%IPrint .GE. 12) then
           call lsheader(lupri,'Generalized inverse of B matrix for delocalized internals')
           call output(Bs_inv,1,N_Cart,1,N_Cart-6,N_Cart,N_Cart-6,1,LUPRI)
        Endif
        !
     Endif ! Finished
  Endif ! Converged
Enddo
! Print summary
If (optinfo%IPrint .GE. 5) then
  call lsheader(lupri,'Back transformation converged in:')
  WRITE(LUPRI,'(I10)') Counter+Cycles
  call lsheader(lupri,'runs.')
  !
  call lsheader(lupri,'Back transformation converged in:')
  WRITE(LUPRI,'(I10)') Cycles+1
  call lsheader(lupri,'cycles.')
Endif
! Copy coordinates back
optinfo%Coordinates = Ini_coord
If (.NOT. optinfo%IBT) then
  ! Get final internals
  Call Get_internals(N_int,N_Cart,CSTEP(1:N_Cart),q,lupri,optinfo)
  Cart_step = D0
  ! Get old internals
  Call Get_internals(N_int,N_Cart,Cart_step,q0,lupri,optinfo)
  ! Printout
  ! Print CSTEP
  call lsheader(lupri,'Cartesian step from back-transformation')
  If (optinfo%oldIBT) call lsheader(lupri,'simplified IBT')
  If (optinfo%New_stepping) call lsheader(lupri,'HOBIT')
   call output(CSTEP(1:N_Cart),1,1,1,N_Cart,1,N_Cart,1,LUPRI)
   ! Modify optinfo%STPINT and optinfo%stpsym
   Do i= 1, N_int
     optinfo%STPINT(i) = q(i) - q0(i)
     If (optinfo%INTCRD(i,1) .GT. 10) optinfo%STPINT(i) = MOD(optinfo%STPINT(i),2.0E0_realk*PI)
     If ((optinfo%INTCRD(i,1) .GT. 20) .AND. ((ABS(optinfo%STPINT(I))-PI) .GT. 0.0E0_realk)) THEN
        If (optinfo%STPINT(i) .GT. 0.0E0_realk) THEN
           optinfo%STPINT(i) = optinfo%STPINT(i) - 2.0E0_realk*PI
        Else
           optinfo%STPINT(i) = optinfo%STPINT(i) + 2.0E0_realk*PI
        Endif
     Endif
  Enddo
  !
  optinfo%STPSYM = CSTEP
  ! Print step in internals
  call lsheader(lupri,'Internal step from back-transformation:')
  If (optinfo%oldIBT) call lsheader(lupri,'simplified IBT')
  If (optinfo%New_stepping) call lsheader(lupri,'HOBIT')
  call output(optinfo%STPINT,1,1,1,N_int,1,N_int,1,LUPRI)
  Int_step_norm = SQRT(DDOT(N_Int,optinfo%STPINT,1,optinfo%STPINT,1))
  WRITE(LUPRI,'(/A,F17.12)')' Norm of internal step', Int_step_norm
  !
  ! Transform step from redundant to delocs
  !
  Call DGEMV('T',N_Int,N_Cart-6,1.0E0_realk,Vectors,N_Int,&
     & optinfo%STPINT(1:N_int),1,1.0E0_realk,Del_step,1)
  call lsheader(lupri,'Step in delocalized internals:')
  call output(Del_step,1,1,1,N_Cart-6,1,N_Cart-6,1,LUPRI)
! IBT
Endif
! 
! Deallocate
!
Call mem_dealloc(Ini_coord)
Call mem_dealloc(Cart_step)
Call mem_dealloc(Int_step)
Call mem_dealloc(Del_step)
Call mem_dealloc(Total_ds)
Call mem_dealloc(Left_ds)
Call mem_dealloc(q0)
Call mem_dealloc(q)
Call mem_dealloc(Bs_inv)
Call mem_dealloc(Vectors)
Call mem_dealloc(B_mat)
Call mem_dealloc(Bs)
! 
End subroutine Back_transform
!================!
! First_inverse  !
!================!
! Calculates vectors and the first inverse
Subroutine First_inverse(optinfo,Bs_inv,vectors,N_Cart,N_Int,lupri)
Implicit none
Integer :: lupri,N_Cart,N_Int
Type(opt_setting) :: optinfo
Real(realk) :: Vectors(N_Int,N_Cart-6)
Real(realk) :: Bs_Inv(N_Cart,N_Cart-6)
Real(realk), pointer :: Left(:,:), Right(:,:),Sing_val(:),Work(:),B_mat(:,:),Sigma(:,:)
Real(realk) :: TS,TE
Integer :: LWork,i,j,Info,N_deloc
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
! Initialize
LWork = 5*N_Int*N_Int
! Allocate singular values, right and left singular vectors, and work
Call mem_alloc(Sing_val,min(N_Int,N_Cart))
Call mem_alloc(Left,N_Int,N_Int)
Call mem_alloc(Right,N_Cart,N_Cart)
Call mem_alloc(Work,LWork)
! Allocate B 
Call mem_alloc(B_mat,N_Int,N_Cart)
! Initialize them
Work = D0 
Left = D0
Right = D0
Sing_val = D0
B_mat = D0

call LSTIMER('START ',TS,TE,lupri)

! Get the B matrix
Call B_matrix(N_int,N_cart,2,B_mat,optinfo)
If (optinfo%IPrint .GE. 12) then
     call lsheader(lupri,'Wilson B matrix [B(ij) = dq(i)/dx(j)] from C++')
     call output(B_mat,1,N_Int,1,N_Cart,N_Int,N_Cart,1,LUPRI)
Endif
! Call LAPACK to do the singular value decomposition
Call DGESVD('A','A',N_Int,N_Cart,B_mat,N_Int, &
& Sing_val,Left,N_Int,Right,N_Cart,Work,LWork,Info)
! Print singular vectors
If (optinfo%IPrint .GE. 12) then 
     call lsheader(lupri,'Singular values')
     Do i = 1, N_Cart
        Write(LUPRI,'(F20.15)') Sing_val(i)
     Enddo
Endif
!
!  Sort vectors to find N_deloc
!
N_deloc = 0
Do i=1, min(N_Cart,N_Int)
   If (abs(Sing_val(i)) .GT. 0.0000000001E0_realk) N_Deloc=N_Deloc+1 
Enddo
Write(*,*)'N DELOC=',N_deloc
! Check if the coordinates are well defined
If ((N_Cart-6) .NE. N_deloc) call lsquit('Internal coordinates ill defined!',lupri)
Bs_inv = D0
call LSTIMER('SVD    ',TS,TE,lupri)
!
! Get pseudo-inverse of Bs:
! Bs(+)=V*Sigma(+)U(*), Bs=U*Sigma*V(*) - SVD 
! Simplified: Bs(+) = V*Sigma(+)
!

! Allocate Sigma
Call mem_alloc(Sigma,N_deloc,N_cart)
Sigma = D0
!
Do i = 1, N_Deloc
   If (abs(Sing_val(i)) .GT. 0.001E0_realk) &
   & Sigma(i,i) = D1/Sing_val(i) 
Enddo
!
!  Do V*Sigma(+) = Bs(+)
!
Call DGEMM('T','T',N_Cart,N_deloc,N_Cart,1.0E0_realk,&
& Right,N_Cart,Sigma,N_deloc,0.0E0_realk,Bs_inv,N_Cart)

! Deallocate Sigma
Call mem_dealloc(Sigma)
! Some printout
If (optinfo%IPrint .GE. 12) then
   call lsheader(lupri,'Generalized inverse of B matrix for delocalized internals')
   call output(Bs_inv,1,N_Cart,1,N_Deloc,N_Cart,N_Deloc,1,LUPRI)
Endif
!
! Copy Left to Vectors
! 
Vectors = Left(:,1:N_deloc)
! Deallocate memory
Call mem_dealloc(Sing_val)
Call mem_dealloc(Left)
Call mem_dealloc(Right)
Call mem_dealloc(Work)
Call mem_dealloc(B_mat)
! 
End subroutine First_inverse

!================!
! Poly_stepping  !
!================!
! Performes a cycle of back transformation
!
Subroutine Poly_stepping(Bs_inv,Vectors,CSTEP,Del_step,N_Cart,N_Int,N_deloc,Conv_thresh,Converged,step_len,lupri,optinfo)
!
Implicit none
Integer :: lupri,N_Cart,N_Int
Type(opt_setting) :: optinfo
Real(realk) :: Vectors(N_Int,N_Deloc)
Real(realk) :: Bs_Inv(N_Cart,N_Deloc)
Real(realk) :: Del_step(N_Deloc)
Real(realk) :: CSTEP(N_Cart)  ! Cartesian step
Real(realk), pointer :: B_mat(:,:),B(:,:), Int_step(:)
Real(realk) :: NonRed_rms,Cart_rms,Red_rms ! Norm of the residuals
Real(realk), pointer ::  Inv_test(:,:), Sigma(:,:),dx(:), & 
& dq(:),ds(:),Cart_steps(:,:),Residuals(:),Del_residuals(:),q0(:),q(:) 
Real(realk), pointer :: A0(:),A1(:),A2(:),Prev_S_CSTEP(:),S_CSTEP(:), x0(:)
Real(realk) :: TE,TS ! CPU time
Real(realk) :: Conv_thresh, step_len
Integer :: i,j,N_deloc
Logical :: Converged
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
! Common allocations
Call mem_alloc(dx,N_Cart)
Call mem_alloc(dq,N_Int)
Call mem_alloc(q,N_Int)
Call mem_alloc(Int_step,N_Int)
Call mem_alloc(q0,N_Int)
Call mem_alloc(Residuals,N_Int)
Call mem_alloc(Del_residuals,N_Deloc)
Call mem_alloc(ds,N_Deloc)
Call mem_alloc(Cart_steps,optinfo%Deriv_order,N_cart)
! For Shanks transformation
If (optinfo%Shanks) then
   Call mem_alloc(A0,N_cart)
   Call mem_alloc(A1,N_cart)
   Call mem_alloc(A2,N_cart)
   Call mem_alloc(Prev_S_CSTEP,N_cart)
   Call mem_alloc(S_CSTEP,N_cart)
   Call mem_alloc(x0,N_cart)
! Get x0
   Do i = 1,N_Cart/3
      x0(3*i-2:3*i) = optinfo%Coordinates(:,i)
   Enddo
Endif
!
Cart_steps = D0
! We find the step in primitive internals corresponding to Del_step
Call DGEMV('N',N_Int,N_deloc,1.0E0_realk,Vectors,N_Int,Del_step,1,&
& 0.0E0_realk, Int_step,1)
!First estimate of Cartesian step
dx = D0
ds = D0
q0 = D0
CSTEP = D0
Call DGEMV('N',N_cart,N_deloc,1.0E0_realk,Bs_inv,N_Cart,Del_step,1,&
& 0.0E0_realk,dx,1)
Cart_steps(1,:) = dx !*step_len
CSTEP = dx !*step_len
dx = D0
! Shanks
If (optinfo%Shanks) then
   A0 = x0
   A1 = A0 + CSTEP
   S_CSTEP = A1
Endif
! Get old internals
Call Get_internals(N_int,N_Cart,dx,q0,lupri,optinfo)

! Calculate the residuals
Call Get_residuals(N_int,N_Cart,optinfo,Residuals,q0,Int_step,CSTEP,q,lupri)
Call DGEMV('T',N_int,N_deloc,1.0E0_realk,Vectors,N_int,&
   & Residuals,1,0.0E0_realk,Del_residuals,1)
! Calculate RMS
NonRed_rms =(D1/SQRT(Real(N_deloc)))*SQRT(DDOT(N_deloc,Del_residuals,1,Del_residuals,1))
Red_rms = (D1/SQRT(Real(N_int)))*SQRT(DDOT(N_Int,Residuals,1,Residuals,1))
Cart_rms = (D1/SQRT(Real(N_Cart)))*SQRT(DDOT(N_Cart,Cart_steps(1,:),1,Cart_steps(1,:),1))
! Printout
If (optinfo%IPrint .GE. 5) then
   call lsheader(lupri,'Residuals for nonredundant coordinates,order')
   WRITE(LUPRI,'(A,I10)') 'Order:    ', 1
   call output(Del_residuals,1,1,1,N_Deloc,1,N_Deloc,1,LUPRI)
   WRITE(LUPRI,'(/A,F20.15)')' RMS of the residuals', NonRed_rms
   call lsheader(lupri,'Residuals for redundant coordinates,order')
   WRITE(LUPRI,'(A,I10)') 'Order:    ', 1
   call output(Residuals,1,1,1,N_int,1,N_int,1,LUPRI)
   WRITE(LUPRI,'(/A,F20.15)')' RMS of the redundant residuals', Red_rms
   call lsheader(lupri,'Residuals for Cartesian coordinates,order')
   WRITE(LUPRI,'(A,I10)') 'Order:    ', 1
   call output(Cart_steps(1,:),1,1,1,N_Cart,1,N_Cart,1,LUPRI)
   WRITE(LUPRI,'(/A,F17.12)')' RMS of the Cartesian residuals', Cart_rms
   call lsheader(lupri,'Cartesian step from high-order derivatives,first estimate')
   call output(CSTEP,1,1,1,N_Cart,1,N_Cart,1,LUPRI)
Endif
! A loop over orders of differentiation
Do i =2, optinfo%Deriv_order
   If (optinfo%Shanks) Prev_S_CSTEP = S_CSTEP
   ! Find dq
   call LSTIMER('START ',TS,TE,lupri)
   Call Int_deriv(optinfo,i,N_Int,N_Cart,dq,optinfo%Deriv_order,Cart_steps)
   call LSTIMER('Internal der.',TS,TE,lupri)
   ! Transform dq to delocs
   Call DGEMV('T',N_int,N_deloc,1.0E0_realk,Vectors,N_int,&
      & dq,1,0.0E0_realk,ds,1)
   ! Get new Cartesian step correction 
   Call DGEMV('N',N_Cart,N_deloc,-1.0E0_realk,Bs_inv,N_Cart,ds,&
   & 1,0.0E0_realk,dx,1)
   ! Update Cartesian step
   Cart_steps(i,:) = dx !*(step_len**i)
   !
   CSTEP = CSTEP + dx !*(step_len**i)
   ! Do Shanks transformation if requested
   If (optinfo%Shanks) then
      A2 = x0 + CSTEP
      S_CSTEP = Shanks(A0,A1,A2,N_Cart)
      A0 = A1
      A1 = A2
      dx = S_CSTEP - Prev_S_CSTEP
      Write(*,*)'Shanks=',S_CSTEP(1)
      Write(*,*)'Step=  ',CSTEP(1) + x0(1)
   Endif
   ! Calculate the residuals
   If (optinfo%Shanks) then
      Call Get_residuals(N_int,N_Cart,optinfo,Residuals,q0,Int_step,S_CSTEP-x0,q,lupri)
   Else
      Call Get_residuals(N_int,N_Cart,optinfo,Residuals,q0,Int_step,CSTEP,q,lupri)
   Endif
   Call DGEMV('T',N_int,N_deloc,1.0E0_realk,Vectors,N_int,&
    & Residuals,1,0.0E0_realk,Del_residuals,1)
   ! Calculate RMS
   NonRed_rms =(D1/SQRT(Real(N_deloc)))*SQRT(DDOT(N_deloc,Del_residuals,1,Del_residuals,1))
   Red_rms = (D1/SQRT(Real(N_int)))*SQRT(DDOT(N_Int,Residuals,1,Residuals,1))
   Cart_rms = (D1/SQRT(Real(N_Cart)))*SQRT(DDOT(N_Cart,dx,1,dx,1))
   ! Printout
   If (optinfo%IPrint .GE. 5) then
      call lsheader(lupri,'Residuals for nonredundant coordinates,order')
      WRITE(LUPRI,'(A,I10)') 'Order:    ', i
      call output(Del_residuals,1,1,1,N_Deloc,1,N_Deloc,1,LUPRI)
      WRITE(LUPRI,'(/A,F20.15)')' RMS of the residuals', NonRed_rms
      call lsheader(lupri,'Residuals for redundant coordinates,order')
      WRITE(LUPRI,'(A,I10)') 'Order:    ', i
      call output(Residuals,1,1,1,N_int,1,N_int,1,LUPRI)
      WRITE(LUPRI,'(/A,F20.15)')' RMS of the redundant residuals', Red_rms
      WRITE(LUPRI,'(A,I10)') 'Order:    ', i
      call lsheader(lupri,'Residuals for Cartesian coordinates,order')
      WRITE(LUPRI,'(A,I10)') 'Order:    ', i
      call output(dx,1,1,1,N_Cart,1,N_Cart,1,LUPRI)
      WRITE(LUPRI,'(/A,F20.15)')' RMS of the Cartesian residuals', Cart_rms
      call lsheader(lupri,'Cartesian step from high-order derivatives')
      call output(CSTEP,1,1,1,N_Cart,1,N_Cart,1,LUPRI)
   Endif
 ! Check convergence
   If ((NonRed_rms .LE. Conv_thresh) .AND. (Cart_rms .LE. Conv_thresh)) then 
      Converged = .TRUE.
   Endif
 ! Convergence not obtained and maximal derivative order was reached 
   If (Converged) EXIT
 ! 
Enddo
! Print x(k) vectors
If (optinfo%IPrint .GE. 12) then
   Do i = 1, optinfo%Deriv_order
      call lsheader(lupri,'x**k vector')
      WRITE(LUPRI,'(A,I10)') 'Order:    ', i
      Do j = 1, N_Cart
         Write(LUPRI,'(F20.15)') Cart_steps(i,j) 
      Enddo
   Enddo
Endif
! Copy step for Shanks transformation
If (optinfo%Shanks) CSTEP = S_CSTEP - x0
! Deallocate memory
Call mem_dealloc(dx)
Call mem_dealloc(dq)
Call mem_dealloc(q)
Call mem_dealloc(q0)
Call mem_dealloc(ds)
Call mem_dealloc(Int_step)
Call mem_dealloc(Cart_steps)
Call mem_dealloc(Residuals)
Call mem_dealloc(Del_Residuals)
! For Shanks transformation
If (optinfo%Shanks) then
   Call mem_dealloc(A0)
   Call mem_dealloc(A1)
   Call mem_dealloc(A2)
   Call mem_dealloc(Prev_S_CSTEP)
   Call mem_dealloc(S_CSTEP)
   Call mem_dealloc(x0)
Endif
!
End subroutine Poly_stepping
!==========!
! IBT      !
!==========!
! Performs iterative back-transformation
! recalculating B and B(+) each iteration
! It is called after a linear approximation to x is performed
!
Subroutine IBT(Bs_inv,N_Cart,N_Int,N_del,optinfo,Cart_step,Del_step,Del_trans,q0,lupri,Converged,Conv_thresh)
Implicit none
Real(realk) :: Cart_step(N_Cart)
Real(realk) :: q0(N_int), Del_step(N_del) 
Real(realk) :: Bs_inv(N_Cart,N_del)
Real(realk) :: Del_trans(N_int,N_Del)    ! Transformation from redundant to non-redundant
Real(realk), pointer :: B_mat(:,:), Bs(:,:), Ini_coord(:,:),dx(:), Residuals(:), &
& Del_residuals(:),q(:)
Type(opt_setting) :: optinfo
Real(realk), pointer :: Int_step(:)
Logical :: Converged
Real(realk) :: NonRed_rms,Cart_rms,Red_rms, Int_step_norm, Conv_thresh 
Integer :: i, j, lupri, N_del, N_int, N_Cart
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
! Allocate
call mem_alloc(B_mat,N_int,N_Cart)
call mem_alloc(Bs,N_del,N_Cart)
call mem_alloc(q,N_int)
call mem_alloc(Int_step,N_int)
call mem_alloc(Residuals,N_int)
call mem_alloc(Ini_Coord,3,N_Cart)
call mem_alloc(dx,N_Cart)
call mem_alloc(Del_residuals,N_del)
! We find the step in primitive internals corresponding to Del_step
Call DGEMV('N',N_Int,N_del,1.0E0_realk,Del_trans,N_Int,Del_step,1,&
& 0.0E0_realk, Int_step,1)
!First estimate of Cartesian step
dx = D0
q0 = D0
Cart_step = D0
Call DGEMV('N',N_cart,N_del,1.0E0_realk,Bs_inv,N_Cart,Del_step,1,&
& 0.0E0_realk,dx,1)
Cart_step = dx !*step_len
dx = D0
! Get old internals
Call Get_internals(N_int,N_Cart,dx,q0,lupri,optinfo)

! Calculate the residuals
Call Get_residuals(N_int,N_Cart,optinfo,Residuals,q0,Int_step,Cart_step,q,lupri)
Call DGEMV('T',N_int,N_del,1.0E0_realk,Del_trans,N_int,&
   & Residuals,1,0.0E0_realk ,Del_residuals,1)
! Calculate RMS
NonRed_rms =(D1/SQRT(Real(N_del)))*SQRT(DDOT(N_del,Del_residuals,1,Del_residuals,1))
Red_rms = (D1/SQRT(Real(N_int)))*SQRT(DDOT(N_Int,Residuals,1,Residuals,1))
Cart_rms = (D1/SQRT(Real(N_Cart)))*SQRT(DDOT(N_Cart,Cart_step,1,Cart_step,1))
! Printout
If (optinfo%IPrint .GE. 5) then
   call lsheader(lupri,'Residuals for nonredundant coordinates,order')
   WRITE(LUPRI,'(A,I10)') 'Order:    ', 1
   call output(Del_residuals,1,1,1,N_Del,1,N_Del,1,LUPRI)
   WRITE(LUPRI,'(/A,F20.15)')' RMS of the residuals', NonRed_rms
   call lsheader(lupri,'Residuals for redundant coordinates,order')
   WRITE(LUPRI,'(A,I10)') 'Order:    ', 1
   call output(Residuals,1,1,1,N_int,1,N_int,1,LUPRI)
   WRITE(LUPRI,'(/A,F20.15)')' RMS of the redundant residuals', Red_rms
   call lsheader(lupri,'Residuals for Cartesian coordinates,order')
   WRITE(LUPRI,'(A,I10)') 'Order:    ', 1
   call output(Cart_step,1,1,1,N_Cart,1,N_Cart,1,LUPRI)
   WRITE(LUPRI,'(/A,F17.12)')' RMS of the Cartesian residuals', Cart_rms
   call lsheader(lupri,'Cartesian step from high-order derivatives,first estimate')
   call output(Cart_step,1,1,1,N_Cart,1,N_Cart,1,LUPRI)
Endif
!
Residuals = D0
Del_residuals = D0
! Copy step
dx = Cart_step
! Copy initial coordinates
Ini_coord = optinfo%Coordinates
!
! Do the loop for reserved number of times
!
Do j = 2, optinfo%Deriv_order
   ! We get new coordinates and B_matrix
   Do i = 1, N_Cart/3
      optinfo%Coordinates(:,i) = optinfo%Coordinates(:,i) + dx(i*3-2:i*3) 
   Enddo
 If (optinfo%IBT) then ! We recalculate Bs(+)
   !
   Call B_matrix(N_int,N_cart,2,B_mat,optinfo)
   !
   ! We get the Bs matrix for delocalized internals
   Call DGEMM('T','N',N_del,N_Cart,N_int,1.0E0_realk,&
   & Del_trans,N_int,B_mat,N_Int,0.0E0_realk,Bs,N_del)
   ! Get pseudoinverse for the B_s
   Call Get_B_Inv(Bs,Bs_inv,N_del,N_Cart)
 Endif ! IBT
   ! 
   !  Get new Cartesian displacement and update the total step
   !
   Call DGEMV('N',N_cart,N_del,1.0E0_realk,Bs_inv,N_Cart,Del_residuals,1,&
   & 0.0E0_realk,dx,1)
   Cart_step = Cart_step + dx
   ! Calculate the residuals
   Call Get_residuals(N_int,N_Cart,optinfo,Residuals,q0,Int_step,dx,q,lupri)
   Call DGEMV('T',N_int,N_del,1.0E0_realk,Del_trans,N_int,&
   & Residuals,1,0.0E0_realk,Del_residuals,1)
   ! Calculate RMS
   NonRed_rms =(D1/SQRT(Real(N_del)))*SQRT(DDOT(N_del,Del_residuals,1,Del_residuals,1))
   Red_rms = (D1/SQRT(Real(N_int)))*SQRT(DDOT(N_Int,Residuals,1,Residuals,1))
   Cart_rms = (D1/SQRT(Real(N_Cart)))*SQRT(DDOT(N_Cart,dx,1,dx,1))
   !
   ! Printout
   If (optinfo%IPrint .GE. 5) then
      call lsheader(lupri,'Residuals for nonredundant coordinates,iteration:')
      WRITE(LUPRI,'(A,I10)') 'Order:    ', j
      call output(Del_residuals,1,1,1,N_Del,1,N_Del,1,LUPRI)
      WRITE(LUPRI,'(/A,F25.18)')'RMS of the residuals', NonRed_rms
      call lsheader(lupri,'Residuals for redundant coordinates,order')
      WRITE(LUPRI,'(A,I10)') 'Order:    ', j
      call output(Residuals,1,1,1,N_int,1,N_int,1,LUPRI)
      WRITE(LUPRI,'(/A,F25.18)')' Norm the residuals', Red_rms
      call lsheader(lupri,'Residuals for Cartesian coordinates,order')
      WRITE(LUPRI,'(A,I10)') 'Order:    ', j
      call output(dx,1,1,1,N_Cart,1,N_Cart,1,LUPRI)
      WRITE(LUPRI,'(/A,F25.18)')' RMS of the Cartesian residuals', Cart_rms
      call lsheader(lupri,'Cartesian step from a proper IBT, not finished')
      call output(Cart_step,1,1,1,N_Cart,1,N_Cart,1,LUPRI)
   Endif
   ! Check convergence
   If ((NonRed_rms .LE. Conv_thresh) .AND. (Cart_rms .LE. Conv_thresh)) then 
      Converged = .TRUE.
   Endif
   !
   If (Converged) EXIT
   ! 
Enddo
! Final printout and modifications
! Print Cart_step
call lsheader(lupri,'Cartesian step from iterative back transformation')
call output(Cart_step(1:N_Cart),1,1,1,N_Cart,1,N_Cart,1,LUPRI)
! Modify optinfo%STPINT and optinfo%stpsym
Do i= 1, N_int
    optinfo%STPINT(i) = q(i) - q0(i)
    If (optinfo%INTCRD(i,1) .GT. 10) optinfo%STPINT(i) = MOD(optinfo%STPINT(i),2.0E0_realk*PI)
    If ((optinfo%INTCRD(i,1) .GT. 20) .AND. ((ABS(optinfo%STPINT(I))-PI) .GT. 0.0E0_realk)) THEN
       If (optinfo%STPINT(i) .GT. 0.0E0_realk) THEN
          optinfo%STPINT(i) = optinfo%STPINT(i) - 2.0E0_realk*PI
       Else
          optinfo%STPINT(i) = optinfo%STPINT(i) + 2.0E0_realk*PI
       Endif
    Endif
Enddo
!
optinfo%STPSYM = Cart_step
! Print step in intertnals
call lsheader(lupri,'Internal step from iterative back transformation')
call output(optinfo%STPINT,1,1,1,N_int,1,N_int,1,LUPRI)
Int_step_norm = SQRT(DDOT(N_Int,optinfo%STPINT,1,optinfo%STPINT,1))
WRITE(LUPRI,'(/A,F17.12)')' Norm of internal step', Int_step_norm
!
! Copy coordinates back
optinfo%Coordinates = Ini_coord
!
! Deallocate everything
!
call mem_dealloc(B_mat)
call mem_dealloc(Bs)
call mem_dealloc(q)
call mem_dealloc(Int_step)
call mem_dealloc(Residuals)
call mem_dealloc(Ini_Coord)
call mem_dealloc(dx)
call mem_dealloc(Del_residuals)
!
End subroutine IBT
!===========!
! Get_B_inv !
!===========!
Subroutine Get_B_inv(B,B_Inv,M,N)
! Calculates generalized inverse of B(M,N) via SVD
Implicit none
Real(realk) :: B(M,N)
Real(realk) :: B_Inv(N,M)
Real(realk), pointer :: Right(:,:),Left(:,:),Sing_val(:),Work(:), &
Inv_test(:,:), Sigma(:,:)
Integer :: M,N,Lwork,i,Info
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
!
LWork = 5*M*N
! Allocate memory
Call mem_alloc(Right,N,N)
Call mem_alloc(Left,M,M)
Call mem_alloc(Work,LWork)
Call mem_alloc(Sing_val,min(N,M))
Call mem_alloc(Inv_test,N,M)
Call mem_alloc(Sigma,M,N)
Sigma = D0
! Call SVD from LAPACK
Inv_test = D0
Call DGESVD('A','A',M,N,B,M, &
& Sing_val,Left,M,Right,N,Work,LWork,Info)
! Allocate Sigma
Sigma = D0
!
Do i = 1, min(M,N)
   If (abs(Sing_val(i)) .GT. 0.001E0_realk) &
      & Sigma(i,i) = D1/Sing_val(i) 
Enddo
Call DGEMM('T','T',N,M,N,1.0E0_realk,&
& Right,N,Sigma,M,0.0E0_realk,Inv_test,N)
Call DGEMM('N','T',N,M,M,1.0E0_realk,&
& Inv_test,N,Left,M,0.0E0_realk,B_Inv,N)
! Deallocate memory
Call mem_dealloc(Right)
Call mem_dealloc(Left)
Call mem_dealloc(Work)
Call mem_dealloc(Sing_val)
Call mem_dealloc(Inv_test)
Call mem_dealloc(Sigma)
!
End subroutine Get_B_inv
!
!=================!
! Get_internals   !
!=================!
! Call cpp program to calculate internals
! at modified geometry
Subroutine Get_internals(N_internal,N_cartesian,Cart_step,New_internals,lupri,optinfo)
Implicit none
Type(opt_setting) :: optinfo
INTEGER :: ldx,i_type,N_internal,N_cartesian,i,j,k,lupri
Real(realk), dimension(N_Cartesian) :: Cart_step
Real(realk), dimension(N_Internal) :: New_internals
Real(realk), pointer :: Cartesian(:,:)
Real(realk), pointer :: q_Internal(:)
Integer, dimension(4) :: i_centre
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
!
New_internals = D0
ldx = 2
Call mem_alloc(Cartesian,2,12)
Call mem_alloc(q_Internal,2)
!
Do i=1, optinfo%NIntCoord
   Cartesian = D0
   q_internal = D0
   ! First, bond lengths
   If (optinfo%INTCRD(i,1) .LT. 10) i_type=1
   ! Second, angles
   If ((optinfo%INTCRD(i,1) .GT. 10) .AND. (optinfo%INTCRD(i,1) .LT. 20)) i_type=2
   ! Third, dihedrals
   If (optinfo%INTCRD(i,1) .GT. 20) i_type=3
   ! Get Cartesians
   Do j = 1,(i_type+1)
      i_centre(j) = optinfo%INTCRD(i,j+1)
      Cartesian(1,(j*3-2):(j*3)) = optinfo%Coordinates(:,i_centre(j)) + &
      & Cart_step((i_centre(j)*3-2):i_centre(j)*3)
   Enddo
   ! Get the internals
   call dqdx(i_type,1,q_Internal,Cartesian,ldx)   
   !  Fill in the new internals
   New_internals(i) = q_Internal(1)
   !
Enddo
! Printout
If (optinfo%IPrint .GE. 10) then
   call lsheader(lupri,'Internals from Get_internals')
   call output(New_internals,1,1,1,N_internal,1,N_internal,1,LUPRI)
Endif
! Deallocate memory
Call mem_dealloc(Cartesian)
Call mem_dealloc(q_Internal)

End subroutine Get_internals
!
!==========!
! B_matrix !
!==========!
! Calls cpp program to get B matrix
Subroutine B_matrix(N_internal,N_cartesian,i_order,B,optinfo)
Implicit none
Type(opt_setting) :: optinfo
INTEGER :: ldx,i_order,max_order,i_type,N_internal,N_cartesian,i,j,k
Real(realk) :: B(N_internal,N_cartesian),B_row(12) ! Wilson's B matrix
Real(realk), pointer :: Cartesian(:,:)
Real(realk), pointer :: q_Internal(:)
Integer, dimension(4) :: i_centre
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
!
ldx = i_order + 1
Call mem_alloc(Cartesian,ldx,12)
Call mem_alloc(q_Internal,i_order+1)
B_row = D0
B = D0
!
Do i=1, optinfo%NIntCoord
   Cartesian = 0.0E0_realk
   q_internal = 0.0E0_realk
   ! First, bond lengths
   If (optinfo%INTCRD(i,1) .LT. 10) i_type=1
   ! Second, angles
   If ((optinfo%INTCRD(i,1) .GT. 10) .AND. (optinfo%INTCRD(i,1) .LT. 20)) i_type=2
   ! Third, dihedrals
   If (optinfo%INTCRD(i,1) .GT. 20) i_type=3
   ! Get Cartesians
   Do j = 1, i_type+1
      i_centre(j) = optinfo%INTCRD(i,j+1)
      Cartesian(1,(j*3-2):(j*3)) = optinfo%Coordinates(:,i_centre(j)) 
   Enddo
   ! Get the derivatives
   Do j = 1, (i_type+1)*3
      Cartesian(2,j) = Cartesian(2,j) + 1.0E0_realk
      call dqdx(i_type,i_order,q_Internal,Cartesian,ldx)   
      Cartesian(2,j) = Cartesian(2,j) - 1.0E0_realk
      B_row(j) = q_internal(2)
   Enddo
   !  Fill in the B matrix
   Do j = 1,(i_type+1)
      B(i,(i_centre(j)*3-2):(i_centre(j)*3)) = B_row((j*3-2):(j*3))
   Enddo
   !
Enddo
! Deallocate memory
Call mem_dealloc(Cartesian)
Call mem_dealloc(q_Internal)
!
End subroutine B_matrix
!=============!
! Int_deriv   !
!=============!
! Calls cpp program to get dq/dx
Subroutine Int_deriv(optinfo,i_order,N_internal,N_Cartesian,dq,der_ord,Cart_steps)
Implicit none
Type(opt_setting) :: optinfo
INTEGER :: ldx,i_order,i_type,N_internal,N_Cartesian,i,j,k,der_ord
Real(realk), pointer :: q_Internal(:)
Real(realk), dimension(N_internal) :: dq
Real(realk), pointer :: Cartesian(:,:)
Real(realk), dimension(der_ord,N_Cartesian) :: Cart_steps ! Contains the
! Cartesian step corrections for all orders
Real(realk), dimension(i_order) :: numerical_errors ! monitor how
!well we solve the inversion problem
Integer, dimension(4) :: i_centre
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
!
ldx = i_order + 1
dq = D0
! Mem_alloc is a Dalton wrapper for allocate, use any alternative
Call mem_alloc(q_Internal,i_order+1)
Call mem_alloc(Cartesian,ldx,12)
!
Do i=1, i_order
   numerical_errors(i) = 0
enddo
!
Do i=1, optinfo%NIntCoord
   Cartesian = D0
   q_internal = D0
   ! optinfo%INTCRD contains information about the types of internals
   ! Instead of using it define yourself i_type.

   ! First, bond lengths
   If (optinfo%INTCRD(i,1) .LT. 10) i_type=1
   ! Second, angles
   If ((optinfo%INTCRD(i,1) .GT. 10) .AND. (optinfo%INTCRD(i,1) .LT. 20)) i_type=2
   ! Third, dihedrals
   If (optinfo%INTCRD(i,1) .GT. 20) i_type=3
   ! Get Cartesians
   Do j = 1,(i_type+1)
      ! Again, define centre(i) - number of atom which forms the internal - as
      ! you like
      i_centre(j) = optinfo%INTCRD(i,j+1)
      ! optinfo%Coordinates is a dalton structure: 3*N_atoms, use any
      ! alternative
      Cartesian(1,(j*3-2):(j*3)) = optinfo%Coordinates(:,i_centre(j))
      ! And Cartesian steps
      Do k = 1,(i_order-1)
         Cartesian(k+1,(j*3-2):(j*3)) = Cart_steps(k,(i_centre(j)*3-2):(i_centre(j)*3))
      Enddo
      !
   Enddo
   ! Get the derivatives: CPP is called
   call dqdx(i_type,i_order,q_Internal,Cartesian,ldx)   
   ! Fill in dq 
   dq(i) = q_Internal(i_order+1)
   !
   do j=3,i_order
      numerical_errors(j) = numerical_errors(j) + q_Internal(j)**2
   enddo
   !
Enddo
!
!do j=3,i_order
!   numerical_errors(j) = sqrt(numerical_errors(j)/optinfo%NIntCoord)
! RMS error
!   print *,'RMS error at order',j-1,' out of ',i_order,'is',numerical_errors(j)
!enddo
! Deallocate memory
Call mem_dealloc(Cartesian)
Call mem_dealloc(q_Internal)
!
End subroutine Int_deriv
!=================!
! Get_Residuals   !
!=================!
Subroutine Get_Residuals(N_internal,N_Cartesian,optinfo,Residuals,q0,dq,Cart_step,q,lupri)
! Computes differences between desired step 
! and updated step in internals
Implicit none
Integer :: i, N_internal,N_Cartesian,lupri
Real(realk), dimension(N_internal) :: q0,q,dq,Residuals 
Real(realk), dimension(N_Cartesian) :: Cart_step 
Real(realk) :: R_norm
Type(opt_setting) :: optinfo
Real, parameter:: D0 = 0.0E0_realk, D1 = 1.0E0_realk
! Get internals from modified Cartesians 
Call Get_internals(N_internal,N_Cartesian,Cart_step,q,lupri,optinfo)
! Printout
If (optinfo%IPrint .GE. 5) then
   call lsheader(lupri,'Internals after approximate step')
   call output(q,1,1,1,N_internal,1,N_internal,1,LUPRI)
Endif
! Estimate residuals
Residuals = dq - (q - q0)
! Correct residuals for bond angles and dihedrals
Do i = 1, N_internal
   !  All angles
   If ((optinfo%INTCRD(i,1) .GT. 10) ) Residuals(i) = MOD(Residuals(i),2.0E0_realk*PI)

!      If ( (q0(i)+dq(i)) .LT. 0.00000000D+0 ) then 
!         Residuals(i) = q(i) + (q0(i) + dq(i))
!      Else
!         If ( (q0(i)+dq(i)) .GT. PI ) then
!            Residuals(i) = q(i) - (2*PI - q0(i) - dq(i))
!         Else
!            Residuals(i) = q(i) - (q0(i) + dq(i))
!         Endif
!      Endif
      !

   ! Dihedrals only
   If ((optinfo%INTCRD(i,1) .GT. 20) .AND. (ABS(Residuals(i)) .GT. PI)) then
       If (Residuals(i) .GT. -0.0001E0_realk) then
          Residuals(i) = Residuals(i) - 2.0E0_realk*PI
       Else
          Residuals(i) = Residuals(i) + 2.0E0_realk*PI
       Endif
    Endif
    ! One more correction to the angles residuals: zero everything close to PI
    If ((optinfo%INTCRD(i,1) .GT. 10) .AND. (ABS(ABS(Residuals(i))-PI) .LE. 0.00001E0_realk ) ) &
       & Residuals(i) = D0
    ! Zero everything close to zero
    IF (ABS(Residuals(i)) .LT. 1.0E-14_realk) Residuals(i) = D0

!      If ( ((q0(i)+dq(i)) .LT. -PI) .AND. (q(i) .GT. 0.00000000D+0) ) then 
!          If (i .EQ. 38) write(*,*) 'Here 111'
!          Residuals(i) = q(i) - (q0(i) + dq(i)) - 2*PI
!      Else
!          If ( ((q0(i)+dq(i)) .GT. PI ) .AND. (q(i) .LT. 0.00000000D+0) ) then
!             If (i .EQ. 38) write(*,*) 'Here 222'
!                 Residuals(i) = q(i) - (q0(i) + dq(i)) + 2*PI
!          Else
!             If (i .EQ. 38) write(*,*) 'Here 333'
!             If (i .EQ. 38) write(*,*) q(i),q0(i),dq(i)
!             Residuals(i) = q(i) - (q0(i) + dq(i)) 
!          Endif
!      Endif
!      !
!   Endif
!

Enddo   ! Loop over internals
!
R_norm = sqrt(dot_product(Residuals,Residuals))
!
End subroutine Get_Residuals
!=============!
!  Shanks     !
!=============!
! A function performing 
! Shanks transformation of vectors
Function Shanks(A0,A1,A2,N) result(S_transform)
Implicit none
Integer :: N,i
Real(realk) :: S_transform(N)
Real(realk) :: A0(N)
Real(realk) :: A1(N)
Real(realk) :: A2(N)
!
Do i=1,N
   S_transform(i) = (A2(i)*A0(i) - A1(i)**2)/(A2(i)-2*A1(i)+A0(i))
Enddo
!
End function Shanks
!
End module q_to_x_mod





















