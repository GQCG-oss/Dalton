module unittest_mat_dense
use unittest_util
use matrix_operations

contains

   subroutine MatrixTestDense
   implicit none
      type(matrix) :: A, B, C, D, E, F, G, H, L, M, N, O, P, Q
      integer      :: i, j, pos, ierr, correct, nnz, io
      real(realk)  :: sumElms, tr, tr2, dot1, dot2, norm, val, sum1, sum2
      real(realk),allocatable :: Afull(:,:), eival1(:), eival2(:), Dfull(:,:),Efull(:,:),Ffull(:,:)
      real(realk),allocatable :: eival_r(:), eival_im(:), eival_denom(:), eorb(:), Gfull(:,:), Hfull(:,:)

      call reset_counters

!=================================
!    Unit test - mat_select_type:
!=================================

      !First, reset matrix type (2 is default):
      call mat_select_type(0,6)
      !Then select dense matrices:
      call mat_select_type(2,6)
      !UNITTEST mat_select_type:
      call AssertEquals(2, matrix_type, "Matrix type not correctly chosen in mat_select_type")

!=================================
!    Unit test - mat_init:
!=================================

      call mat_init(A,5,7)
      !Test nrow and ncol:
      call AssertEquals(5, A%nrow, "Nrow not correctly set in mat_init")
      call AssertEquals(7, A%ncol, "Ncol not correctly set in mat_init")
      !Test if elms pointer is associated:
      call AssertEquals(ASSOCIATED(A%elms),.true.,"Elms pointer not associated in mat_dense_init")
      !Test allocated size of A:
      call AssertEquals(35, size(A%elms), "Allocated size not correct in mat_init")
      
!=================================
!    Unit test - mat_zero:
!=================================

      call mat_zero(A)
      sumElms = 0.0d0
      do i = 1, 35
         sumElms = sumElms + A%elms(i)
      enddo
      call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_dense_zero not working properly")

!======================================
!    Unit test - mat_create_elm:
!======================================

      call mat_create_elm(3,6,1.5d0,A)
      call AssertEquals(1.5d0,A%elms(28),1.0d-20,"Mat_create_elm not working properly")

      call mat_create_elm(1,6,3.4d0,A)
      call mat_create_elm(2,2,5.2d0,A)
      call mat_create_elm(5,4,2.3d0,A)
      call mat_create_elm(5,7,1.4d0,A)

!======================================
!    Unit test - mat_to_full
!======================================

     allocate(Afull(5,7))
     call mat_to_full(A,2.0d0,Afull)

     sumElms = 0.0d0
     do i = 1, 5
        do j = 1, 7
           sumElms = sumElms + Afull(i,j) - 2.0d0*a%elms(a%nrow*(j-1)+i)
           !print *, 'Afull(i,j)', Afull(i,j)
           !print *, 'a%elms(a%nrow*(j-1)+i', a%elms(a%nrow*(j-1)+i)
        enddo
     enddo
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_to_full not working properly")

!======================================
!    Unit test - mat_set_from_full
!======================================

     call mat_zero(A)
     call mat_set_from_full(Afull,0.5d0, A)

     sumElms = 0.0d0
     do i = 1, 5
        do j = 1, 7
           sumElms = sumElms + Afull(i,j) - 2.0d0*a%elms(a%nrow*(j-1)+i)
           !print *, 'Afull(i,j)', Afull(i,j)
           !print *, 'a%elms(a%nrow*(j-1)+i', a%elms(a%nrow*(j-1)+i)
        enddo
     enddo
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_set_from_full not working properly")

     deallocate(Afull)
!======================================
!    Unit test - mat_trans:
!======================================

     call mat_init(B, 7, 5)
     call mat_trans(A,B)
     sumElms = 0.0d0
     do i = 1, 5
        do j = 1, 7
           sumElms = sumElms + b%elms(b%nrow*(i-1)+j) - a%elms(a%nrow*(j-1)+i)
        enddo
     enddo
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_trans not working properly")

!======================================
!    Unit test - mat_assign:
!======================================

     call mat_init(C,5,7)
     call mat_assign(C, A)

     call AssertEquals(A,C,"Mat_assign not working properly")

!======================================
!    Unit test - mat_copy:
!======================================

     call mat_zero(C)
     call mat_copy(2.0d0,A,C)
     sumElms = 0.0d0
     do i = 1, 35
        sumElms = sumElms + 2.0d0*A%elms(i) - C%elms(i)
     enddo
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_copy not working properly")

!======================================
!    Unit test - mat_tr:
!======================================

     call mat_init(D,7,7)

     call mat_create_elm(3,6,1.5d0,D)
     call mat_create_elm(1,6,3.4d0,D)
     call mat_create_elm(2,2,5.2d0,D)
     call mat_create_elm(5,4,2.3d0,D)
     call mat_create_elm(5,7,1.4d0,D)
     call mat_create_elm(7,7,2.4d0,D)

     tr = mat_tr(D)

     sumElms = 0.0d0
     do i = 1, 7
        sumElms = sumElms + d%elms(d%nrow*(i-1)+i) 
        !print *, 'element:', d%nrow*(i-1)+i
        !print *, 'sumElms', sumElms
     enddo

     call AssertEquals(sumElms,tr,1.0d-20,"Mat_tr not working properly")

!======================================
!    Unit test - mat_mul:
!======================================

     !1. Test mat_mul - square matrix - transa = 'n', transb = 'n'

     call mat_init(E,7,7)
     call mat_init(F,7,7)
     call mat_init(G,7,7)

     call mat_create_elm(4,7,2.4d0,E)
     call mat_create_elm(2,7,4.3d0,E)
     call mat_create_elm(3,3,6.1d0,E)
     call mat_create_elm(6,5,3.2d0,E)
     call mat_create_elm(6,1,2.3d0,E)
     call mat_create_elm(1,1,3.3d0,E)

     call mat_mul(D,E,'n','n',4.0d0,0.0d0,F)

     call my_mat_mul(D,E,G)
     !call mat_print(G, 1, 7, 1, 7, 6)

     sumElms = 0.0d0
     do i = 1, 49
        sumElms = sumElms + F%elms(i) - 4.0d0*G%elms(i)
        !print *, 'element:', d%nrow*(i-1)+i
        !print *, 'sumElms', sumElms
     enddo
     
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_mul - square - 'n'/'n' not working properly")

     !2. Test mat_mul - square matrix - transa = 't', transb = 'n'

     call mat_mul(D,E,'t','n',4.0d0,0.0d0,F)

     call mat_init(H,7,7)
     call mat_trans(D,H)
     call my_mat_mul(H,E,G)

     sumElms = 0.0d0
     do i = 1, 49
        sumElms = sumElms + F%elms(i) - 4.0d0*G%elms(i)
     enddo
     
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_mul - square - 't'/'n' not working properly")

     !3. Test mat_mul - square matrix - transa = 'n', transb = 't'

     call mat_mul(D,E,'n','t',4.0d0,0.0d0,F)

     call mat_init(L,7,7)
     call mat_trans(E,L)
     call my_mat_mul(D,L,G)

     sumElms = 0.0d0
     do i = 1, 49
        sumElms = sumElms + F%elms(i) - 4.0d0*G%elms(i)
     enddo
     
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_mul - square - 'n'/'t' not working properly")

     !4. Test mat_mul - square matrix - transa = 't', transb = 't'

     call mat_mul(D,E,'t','t',4.0d0,0.0d0,F)

     call my_mat_mul(H,L,G)

     sumElms = 0.0d0
     do i = 1, 49
        sumElms = sumElms + F%elms(i) - 4.0d0*G%elms(i)
     enddo
     
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_mul - square - 't'/'t' not working properly")

     !5. Test mat_mul - non-square matrix:

     call mat_init(M,5,5)
     call mat_init(N,5,5)

     call mat_mul(A,B,'n','n',6.0d0,0.0d0,M)

     call my_mat_mul(A,B,N)

     sumElms = 0.0d0
     do i = 1, 25
        sumElms = sumElms + M%elms(i) - 6.0d0*N%elms(i)
     enddo
     
     call AssertEquals(sumElms,0.0d0,1.0d-20,"Mat_mul - non-square - not working properly")

!======================================
!    Unit test - mat_TrAB:
!======================================

     tr = mat_trAB(E,F)

     call mat_mul(E,F,'n','n',1.0d0,0.0d0,G)
     tr2 = mat_tr(G) 

     call AssertEquals(tr,tr2,1.0d-20,"Mat_TrAB not working properly")

!======================================
!    Unit test - mat_scal:
!======================================

     !call mat_print(G,1,G%nrow,1,G%ncol,6)

     call mat_scal(0.3d0,G)

     sumElms = 0.0d0
     do i = 1, 49
        sumElms = sumElms + G%elms(i)
     enddo

     call AssertEquals(sumElms,165.3024d0,1.0d-12,"Mat_scal not working properly")

!======================================
!    Unit test - mat_scal_dia:
!======================================

     !call mat_print(G,1,G%nrow,1,G%ncol,6)

     call mat_scal_dia(0.2d0,G)

     sumElms = 0.0d0
     do i = 1, 7
        sumElms = sumElms + G%elms((G%Nrow+1)*i-G%Nrow)
     enddo

     call AssertEquals(sumElms,13.968d0,1.0d-12,"Mat_scal_dia not working properly")

!======================================
!    Unit test - mat_add:
!======================================

     call mat_init(O, 5, 7)
     call mat_init(P, 5, 7)
     call mat_create_elm(4,7,2.4d0,C)
     call mat_create_elm(2,7,4.3d0,C)
     call mat_create_elm(3,3,6.1d0,C)
     call mat_create_elm(5,5,3.2d0,C)
     call mat_create_elm(5,1,2.3d0,C)
     call mat_create_elm(1,1,3.3d0,C)
     !call mat_print(A, 1, A%nrow, 1, A%ncol, 6)
     !call mat_print(C, 1, C%nrow, 1, C%ncol, 6)

     call mat_add(2.0d0, A, 3.0d0, C, O)

     call mat_scal(2.0d0, A)
     call mat_scal(3.0d0, C)
     call my_mat_add(A,C,P)

     call AssertEquals(O,P,"Mat_add not working properly")

!======================================
!    Unit test - mat_daxpy:
!======================================

     call mat_init(Q,5,7)
     call mat_copy(1.0d0, A, O)
     call mat_copy(1.0d0, C, P)

     call mat_daxpy(1.5d0, A, C)

     call mat_scal(1.5d0,O)
     call mat_add(1.0d0, O, 1.0d0, P, Q)

     !call mat_print(Q, 1, Q%nrow, 1, Q%ncol, 6)
     !call mat_print(C, 1, C%nrow, 1, C%ncol, 6)

     call AssertEquals(C,Q,"Mat_daxpy not working properly")

!======================================
!    Unit test - mat_dotproduct:
!======================================

     !call mat_print(A, 1, A%nrow, 1, A%ncol, 6)
     !call mat_print(C, 1, C%nrow, 1, C%ncol, 6)
     dot1 = mat_dotproduct(A,C)

     dot2 = 0.0d0
     do i = 1,a%nrow*a%ncol
        dot2 = dot2 + A%elms(i)*C%elms(i)
     enddo

     call AssertEquals(dot1,dot2,1.0d-12,"Mat_dotproduct not working properly")

!======================================
!    Unit test - mat_sqnorm2:
!======================================

     !call mat_print(A, 1, A%nrow, 1, A%ncol, 6)
     dot1 = mat_dotproduct(A,A)
     dot2 = mat_sqnorm2(A)

     call AssertEquals(dot1,dot2,1.0d-12,"Mat_sqnorm2 not working properly")

!======================================
!    Unit test - mat_outdia_sqnorm2:
!======================================

     norm = mat_outdia_sqnorm2(C)

     call AssertEquals(norm,1971.72d0,1.0d-12,"Mat_outdia_sqnorm2 not working properly")

!======================================
!    Unit test - mat_abs_max_elm:
!======================================

     call mat_create_elm(4,7,-20.4d0,A)
     call mat_create_elm(2,3,-4.3d0,A)
     call mat_create_elm(3,3,-6.1d0,A)
     call mat_create_elm(5,5,-3.2d0,A)
     call mat_create_elm(1,5,-2.3d0,A)
     call mat_create_elm(1,1,-3.3d0,A)
     !call mat_print(A, 1, A%nrow, 1, A%ncol, 6)

     call mat_abs_max_elm(A, val)

     call AssertEquals(val,20.4d0,1.0d-20,"Mat_abs_max_elm not working properly")

!======================================
!    Unit test - mat_max_elm:
!======================================

     call mat_max_elm(A, val)

     call AssertEquals(val,10.4d0,1.0d-20,"Mat_max_elm not working properly")

!======================================
!    Unit test - mat_max_diag_elm:
!======================================

     call mat_create_elm(4,7,-20.4d0,D)
     call mat_create_elm(2,3,-4.3d0,D)
     call mat_create_elm(3,3,-6.1d0,D)
     call mat_create_elm(5,5,-3.2d0,D)
     call mat_create_elm(1,5,-2.3d0,D)
     call mat_create_elm(1,1,-3.3d0,D)
     call mat_max_diag_elm(D, pos, val)

     call AssertEquals(val,-6.1d0,1.0d-20,"Mat_max_diag_elm not working properly (value)")
     call AssertEquals(pos,3,"Mat_max_diag_elm not working properly (position)")

!======================================
!    Unit test - mat_column_norm:
!======================================

     norm = mat_column_norm(D,7,1,7)

     call AssertEquals(norm,423.88d0,1.0d-12,"Mat_column_norm not working properly")

!======================================
!    Unit test - mat_diag_f:
!======================================

     allocate(eival1(7),eival2(7))
     allocate(eival_r(7),eival_im(7),eival_denom(7))
     allocate(Dfull(7,7),Efull(7,7),Ffull(7,7))
       Dfull(1,1)=-32.72931324; Dfull(1,2)=-7.61479236; Dfull(1,3)=-0.01830394; Dfull(1,4)= 0.00000000
       Dfull(2,1)= -7.61479236; Dfull(2,2)=-9.34241438; Dfull(2,3)=-0.21471365; Dfull(2,4)= 0.00000000
       Dfull(3,1)= -0.01830394; Dfull(3,2)=-0.21471365; Dfull(3,3)=-7.54958660; Dfull(3,4)= 0.00000000
       Dfull(4,1)=  0.00000000; Dfull(4,2)= 0.00000000; Dfull(4,3)= 0.00000000; Dfull(4,4)=-7.63234198
       Dfull(6,1)= -1.78002693; Dfull(6,2)=-3.78384114; Dfull(6,3)=-1.57160276; Dfull(6,4)=-2.10798753
       Dfull(7,1)= -1.78002693; Dfull(7,2)=-3.78384114; Dfull(7,3)=-1.57160276; Dfull(7,4)= 2.10798753

       Dfull(1,5)= 0.00000000;  Dfull(1,6)=-1.78002693; Dfull(1,7)=-1.78002693
       Dfull(2,5)= 0.00000000;  Dfull(2,6)=-3.78384114; Dfull(2,7)=-3.78384114
       Dfull(3,5)= 0.00000000;  Dfull(3,6)=-1.57160276; Dfull(3,7)=-1.57160276
       Dfull(4,5)= 0.00000000;  Dfull(4,6)=-2.10798753; Dfull(4,7)= 2.10798753
       Dfull(5,5)=-7.46385840;  Dfull(5,6)= 0.00000000; Dfull(5,7)= 0.00000000
       Dfull(6,5)= 0.00000000;  Dfull(6,6)=-5.09713691; Dfull(6,7)=-1.55333274
       Dfull(7,5)= 0.00000000;  Dfull(7,6)=-1.55333274; Dfull(7,7)=-5.09713691

       Efull(1,1)=1.00000000; Efull(1,2)=0.23670394; Efull(1,3)= 0.00000000; Efull(1,4)= 0.00000000
       Efull(2,1)=0.23670394; Efull(2,2)=1.00000000; Efull(2,3)= 0.00000000; Efull(2,4)= 0.00000000
       Efull(3,1)=0.00000000; Efull(3,2)=0.00000000; Efull(3,3)= 1.00000000; Efull(3,4)= 0.00000000
       Efull(4,1)=0.00000000; Efull(4,2)=0.00000000; Efull(4,3)= 0.00000000; Efull(4,4)= 1.00000000
       Efull(6,1)=0.05490744; Efull(6,2)=0.47954385; Efull(6,3)= 0.22993097; Efull(6,4)= 0.32234014
       Efull(7,1)=0.05490744; Efull(7,2)=0.47954385; Efull(7,3)= 0.22993097; Efull(7,4)=-0.32234014
                                                  
       Efull(1,5)=0.00000000; Efull(1,6)=0.05490744; Efull(1,7)= 0.05490744
       Efull(2,5)=0.00000000; Efull(2,6)=0.47954385; Efull(2,7)= 0.47954385
       Efull(3,5)=0.00000000; Efull(3,6)=0.22993097; Efull(3,7)= 0.22993097
       Efull(4,5)=0.00000000; Efull(4,6)=0.32234014; Efull(4,7)=-0.32234014
       Efull(5,5)=1.00000000; Efull(5,6)=0.00000000; Efull(5,7)= 0.00000000
       Efull(6,5)=0.00000000; Efull(6,6)=1.00000000; Efull(6,7)= 0.24003830
       Efull(7,5)=0.00000000; Efull(7,6)=0.24003830; Efull(7,7)= 1.00000000

     call mat_set_from_full(Dfull,1.0d0,D)
     call mat_set_from_full(Efull,1.0d0,E)

     call mat_diag_f(D,E,eival1,F)
     call RGG(7,7,Dfull,Efull,eival_r,eival_im,eival_denom,1,Ffull,IERR)

     do i = 1,7
        eival_r(i) = eival_r(i)/eival_denom(i)
     enddo

     correct = 0
     do i = 1, 7
        do j = 1,7
           val = eival1(i) - eival_r(j)
           if (abs(val) < 1.0d-12) then
              correct = correct + 1
           endif
        enddo
     enddo

     call AssertEquals(correct,7,"Mat_diag_f not working properly - result differs from RGG")

     deallocate(eival1,eival2)
     deallocate(eival_r,eival_im,eival_denom)
     deallocate(Dfull,Efull,Ffull)

!======================================
!    Unit test - mat_section:
!======================================

     call mat_section(D,1,7,1,5,B)

     val = 0.0d0
     do i = 1, 35
        val = val + D%elms(i) - B%elms(i)
     enddo

     call AssertEquals(val,0.0d0,1.0d-20,"Mat_section not working properly")

!======================================
!    Unit test - mat_section2:
!======================================

     call mat_copy(1.0d0,E,G)

     call mat_section2(B,1,7,1,5,G)

     val = 0.0d0
     do i = 1, 35
        val = val + G%elms(i) - B%elms(i)
     enddo

     do i = 36, 49
        val = val + G%elms(i) - E%elms(i)
     enddo

     call AssertEquals(val,0.0d0,1.0d-20,"Mat_section2 not working properly")

!======================================
!    Unit test - mat_mo_precond:
!======================================

     allocate(eorb(14))
     allocate(dfull(7,7),efull(7,7))

     eorb(1)=-20.2435698066132;       eorb(2)=-1.26098191159298;       eorb(3)=-0.610142947028922
     eorb(4)=-0.451254209418688;      eorb(5)=-0.390335488454913;      eorb(6)=0.591068152739257
     eorb(7)=0.724926445601970;       eorb(8)=2.44583875890112;        eorb(9)=6.87822897159230
     eorb(10)=6.864302064989647E-017; eorb(11)=9.742819426592969E-002; eorb(12)=-2.028335329977485E-018
     eorb(13)=2.64803576008244;       eorb(14)=2.64803576008244

     Dfull = 0.0d0

     Dfull(2,3)=-0.00000001
     Dfull(3,2)= 0.00000001;  Dfull(3,4)=-0.00000005
     Dfull(4,3)= 0.00000005
     Dfull(6,3)= 1.68022125 
     Dfull(7,1)= 0.11192190;  Dfull(7,2)= 0.27256554;  Dfull(7,4)=-1.35122547

     Dfull(1,7)=-0.11192190
     Dfull(2,7)=-0.27256554
     Dfull(3,6)=-1.68022125
     Dfull(4,7)= 1.35122547
     Dfull(6,7)=-0.00000003
     Dfull(7,6)= 0.00000003

     call mat_set_from_full(Dfull,1.0d0,D)

     Efull = 0.0d0

     Efull(2,3)=-0.00000001
     Efull(3,2)= 0.00000001;  Efull(3,4)=-0.00000005
     Efull(4,3)= 0.00000005
     Efull(6,3)= 0.69938633
     Efull(7,1)= 0.00266881;  Efull(7,2)= 0.06862490;  Efull(7,4)=-0.57441238 

     Efull(1,7)=-0.00266881
     Efull(2,7)=-0.06862490
     Efull(3,6)=-0.69938633
     Efull(4,7)= 0.57441238
     Efull(6,7)=-0.00000003
     Efull(7,6)= 0.00000003

     call mat_set_from_full(Efull,1.0d0,E)

     call mat_mo_precond(5,0.0d0,eorb,D)

     call AssertEquals(E,D,"Mat_mo_precond not working properly")
   
     deallocate(eorb)
     deallocate(dfull,efull)

     !call mat_print(D, 1, D%nrow, 1, D%ncol, 6)

!======================================
!    Unit test - mat_ao_precond:
!======================================

     allocate(dfull(7,7),efull(7,7),ffull(7,7),gfull(7,7),hfull(7,7))

       !Dfull = FUP
       Dfull(1,1)=-19.98369891; Dfull(1,2)=-0.35775810; Dfull(1,3)=-0.04811232; Dfull(1,4)= 0.00000000 
       Dfull(2,1)= -0.35775810; Dfull(2,2)=-0.87165022; Dfull(2,3)= 0.05386018; Dfull(2,4)= 0.00000000 
       Dfull(3,1)= -0.04811232; Dfull(3,2)= 0.05386018; Dfull(3,3)=-0.31091742; Dfull(3,4)= 0.00000000 
       Dfull(4,1)=  0.00000000; Dfull(4,2)= 0.00000000; Dfull(4,3)= 0.00000000; Dfull(4,4)=-0.37006383
       Dfull(5,1)=  0.00000000; Dfull(5,2)= 0.00000000; Dfull(5,3)= 0.00000000; Dfull(5,4)= 0.00000000
       Dfull(6,1)=  0.56666555; Dfull(6,2)= 0.35842142; Dfull(6,3)=-0.13961809; Dfull(6,4)=-0.19425934
       Dfull(7,1)=  0.56666555; Dfull(7,2)= 0.35842142; Dfull(7,3)=-0.13961809; Dfull(7,4)= 0.19425934 
                                                                              
       Dfull(1,5)= 0.00000000; Dfull(1,6)= 0.56666555; Dfull(1,7)= 0.56666555
       Dfull(2,5)= 0.00000000; Dfull(2,6)= 0.35842142; Dfull(2,7)= 0.35842142
       Dfull(3,5)= 0.00000000; Dfull(3,6)=-0.13961809; Dfull(3,7)=-0.13961809
       Dfull(4,5)= 0.00000000; Dfull(4,6)=-0.19425934; Dfull(4,7)= 0.19425934
       Dfull(5,5)=-0.28319939; Dfull(5,6)= 0.00000000; Dfull(5,7)= 0.00000000
       Dfull(6,5)= 0.00000000; Dfull(6,6)=-0.30427144; Dfull(6,7)=-0.10032452
       Dfull(7,5)= 0.00000000; Dfull(7,6)=-0.10032452; Dfull(7,7)=-0.30427144

       !Efull = FUQ
       Efull(1,1)= 0.00040467; Efull(1,2)= 0.00679303; Efull(1,3)=-0.00695493; Efull(1,4)= 0.00000000
       Efull(2,1)= 0.00679303; Efull(2,2)= 0.11403140; Efull(2,3)=-0.11674914; Efull(2,4)= 0.00000000 
       Efull(3,1)=-0.00695493; Efull(3,2)=-0.11674914; Efull(3,3)= 0.11953165; Efull(3,4)= 0.00000000 
       Efull(4,1)= 0.00000000; Efull(4,2)= 0.00000000; Efull(4,3)= 0.00000000; Efull(4,4)= 0.27949962
       Efull(5,1)= 0.00000000; Efull(5,2)= 0.00000000; Efull(5,3)= 0.00000000; Efull(5,4)= 0.00000000
       Efull(6,1)= 0.00898455; Efull(6,2)= 0.15081951; Efull(6,3)=-0.15441403; Efull(6,4)=-0.26622324
       Efull(7,1)= 0.00898455; Efull(7,2)= 0.15081951; Efull(7,3)=-0.15441403; Efull(7,4)= 0.26622324 
                                                                             
       Efull(1,5)= 0.00000000; Efull(1,6)= 0.00898455; Efull(1,7)= 0.00898455
       Efull(2,5)= 0.00000000; Efull(2,6)= 0.15081951; Efull(2,7)= 0.15081951
       Efull(3,5)= 0.00000000; Efull(3,6)=-0.15441403; Efull(3,7)=-0.15441403
       Efull(4,5)= 0.00000000; Efull(4,6)=-0.26622324; Efull(4,7)= 0.26622324
       Efull(5,5)= 0.00000000; Efull(5,6)= 0.00000000; Efull(5,7)= 0.00000000
       Efull(6,5)= 0.00000000; Efull(6,6)= 0.45305347; Efull(6,7)=-0.05410153
       Efull(7,5)= 0.00000000; Efull(7,6)=-0.05410153; Efull(7,7)= 0.45305347

       !Ffull = DU
       Ffull(1,1)= 0.99936063; Ffull(1,2)=-0.01073285; Ffull(1,3)= 0.01098864; Ffull(1,4)= 0.00000000 
       Ffull(2,1)=-0.01073285; Ffull(2,2)= 0.81983274; Ffull(2,3)= 0.18446123; Ffull(2,4)= 0.00000000 
       Ffull(3,1)= 0.01098864; Ffull(3,2)= 0.18446123; Ffull(3,3)= 0.81114246; Ffull(3,4)= 0.00000000 
       Ffull(4,1)= 0.00000000; Ffull(4,2)= 0.00000000; Ffull(4,3)= 0.00000000; Ffull(4,4)= 0.64469843
       Ffull(5,1)= 0.00000000; Ffull(5,2)= 0.00000000; Ffull(5,3)= 0.00000000; Ffull(5,4)= 0.00000000
       Ffull(6,1)=-0.01419541; Ffull(6,2)=-0.23829171; Ffull(6,3)= 0.24397098; Ffull(6,4)= 0.33842456
       Ffull(7,1)=-0.01419541; Ffull(7,2)=-0.23829171; Ffull(7,3)= 0.24397098; Ffull(7,4)=-0.33842456 
                                                                             
       Ffull(1,5)= 0.00000000; Ffull(1,6)=-0.01419541; Ffull(1,7)=-0.01419541
       Ffull(2,5)= 0.00000000; Ffull(2,6)=-0.23829171; Ffull(2,7)=-0.23829171
       Ffull(3,5)= 0.00000000; Ffull(3,6)= 0.24397098; Ffull(3,7)= 0.24397098
       Ffull(4,5)= 0.00000000; Ffull(4,6)= 0.33842456; Ffull(4,7)=-0.33842456
       Ffull(5,5)= 1.00000000; Ffull(5,6)= 0.00000000; Ffull(5,7)= 0.00000000
       Ffull(6,5)= 0.00000000; Ffull(6,6)= 0.36248287; Ffull(6,7)= 0.00718130
       Ffull(7,5)= 0.00000000; Ffull(7,6)= 0.00718130; Ffull(7,7)= 0.36248287

       !Gfull = X_AO, input
       Gfull(1,1)= 0.00000000; Gfull(1,2)= 0.23616910; Gfull(1,3)=-0.24142416; Gfull(1,4)= 0.00000000 
       Gfull(2,1)=-0.23616910; Gfull(2,2)= 0.00000000; Gfull(2,3)= 0.00627182; Gfull(2,4)= 0.00000000 
       Gfull(3,1)= 0.24142416; Gfull(3,2)=-0.00627182; Gfull(3,3)= 0.00000000; Gfull(3,4)= 0.00000000 
       Gfull(4,1)= 0.00000000; Gfull(4,2)= 0.00000000; Gfull(4,3)= 0.00000000; Gfull(4,4)= 0.00000000
       Gfull(5,1)= 0.00000000; Gfull(5,2)= 0.00000000; Gfull(5,3)= 0.00000000; Gfull(5,4)= 0.00000000
       Gfull(6,1)=-0.31169810; Gfull(6,2)= 0.01111978; Gfull(6,3)=-0.00308960; Gfull(6,4)=-0.03489263
       Gfull(7,1)=-0.31169810; Gfull(7,2)= 0.01111978; Gfull(7,3)=-0.00308960; Gfull(7,4)= 0.03489263 

       Gfull(1,5)= 0.00000000; Gfull(1,6)= 0.31169810; Gfull(1,7)= 0.31169810
       Gfull(2,5)= 0.00000000; Gfull(2,6)=-0.01111978; Gfull(2,7)=-0.01111978
       Gfull(3,5)= 0.00000000; Gfull(3,6)= 0.00308960; Gfull(3,7)= 0.00308960
       Gfull(4,5)= 0.00000000; Gfull(4,6)= 0.03489263; Gfull(4,7)=-0.03489263
       Gfull(5,5)= 0.00000000; Gfull(5,6)= 0.00000000; Gfull(5,7)= 0.00000000
       Gfull(6,5)= 0.00000000; Gfull(6,6)= 0.00000000; Gfull(6,7)= 0.00000000
       Gfull(7,5)= 0.00000000; Gfull(7,6)= 0.00000000; Gfull(7,7)= 0.00000000

       !Hfull = X_AO, output
       Hfull(1,1)= 0.00000000; Hfull(1,2)= 0.01126235; Hfull(1,3)=-0.01182608; Hfull(1,4)= 0.00000000 
       Hfull(2,1)=-0.01126235; Hfull(2,2)= 0.00000000; Hfull(2,3)= 0.00442884; Hfull(2,4)= 0.00000000 
       Hfull(3,1)= 0.01182608; Hfull(3,2)=-0.00442884; Hfull(3,3)= 0.00000000; Hfull(3,4)= 0.00000000 
       Hfull(4,1)= 0.00000000; Hfull(4,2)= 0.00000000; Hfull(4,3)= 0.00000000; Hfull(4,4)= 0.00000000
       Hfull(5,1)= 0.00000000; Hfull(5,2)= 0.00000000; Hfull(5,3)= 0.00000000; Hfull(5,4)= 0.00000000
       Hfull(6,1)=-0.01502780; Hfull(6,2)= 0.00637966; Hfull(6,3)=-0.00260117; Hfull(6,4)=-0.02480128
       Hfull(7,1)=-0.01502780; Hfull(7,2)= 0.00637966; Hfull(7,3)=-0.00260117; Hfull(7,4)= 0.02480128 
                                                                             
       Hfull(1,5)= 0.00000000; Hfull(1,6)= 0.01502780; Hfull(1,7)= 0.01502780
       Hfull(2,5)= 0.00000000; Hfull(2,6)=-0.00637966; Hfull(2,7)=-0.00637966
       Hfull(3,5)= 0.00000000; Hfull(3,6)= 0.00260117; Hfull(3,7)= 0.00260117
       Hfull(4,5)= 0.00000000; Hfull(4,6)= 0.02480128; Hfull(4,7)=-0.02480128
       Hfull(5,5)= 0.00000000; Hfull(5,6)= 0.00000000; Hfull(5,7)= 0.00000000
       Hfull(6,5)= 0.00000000; Hfull(6,6)= 0.00000000; Hfull(6,7)= 0.00000000
       Hfull(7,5)= 0.00000000; Hfull(7,6)= 0.00000000; Hfull(7,7)= 0.00000000

       call mat_set_from_full(Dfull,1.0d0,D)
       call mat_set_from_full(Efull,1.0d0,E)
       call mat_set_from_full(Ffull,1.0d0,F)
       call mat_set_from_full(Gfull,1.0d0,G)
       call mat_set_from_full(Hfull,1.0d0,H)

       call mat_ao_precond(2,0.0d0,D,E,F,G)

       call AssertEquals(G,H,"Mat_ao_precond not working properly")

       deallocate(dfull,efull,ffull,gfull,hfull)

!=================================
!    Unit test - mat_add_identity:
!=================================

       call mat_add_identity(1.5d0, 2.0d0, L, G)

       sumElms = 0.0d0
       do i = 1, 49
          sumElms = sumElms + G%elms(i)
       enddo

       call AssertEquals(sumElms,53.7d0,1.0d-12,"Mat_add_identity not working properly")

!=================================
!    Unit test - mat_identity:
!=================================

       call mat_identity(L)
       sumElms = 0.0d0
       do i = 1, 49
          sumElms = sumElms + L%elms(i)
       enddo

       call AssertEquals(sumElms,7.0d0,1.0d-12,"Mat_identity not working properly")

!=================================
!    Unit test - mat_get_elm:
!=================================

     !call mat_print(A, 1, A%nrow, 1, A%ncol, 6)
     !call mat_print(B, 1, B%nrow, 1, B%ncol, 6)
     !call mat_print(C, 1, C%nrow, 1, C%ncol, 6)

     sumElms = 0.0d0
     call mat_get_elm (A, 1, 6, val) ; sumElms = val
     call mat_get_elm (B, 2, 2, val) ; sumElms = sumElms + val
     call mat_get_elm (C, 5, 7, val) ; sumElms = sumElms + val

     call AssertEquals(sumElms,10.0575851440430d0,1.0d-12,"Mat_get_elm not working properly")

!=================================
!    Unit test - mat_create_block:
!=================================

     !call mat_print(C, 1, C%nrow, 1, C%ncol, 6)
     allocate(Afull(3,2))
     Afull = 2.35d0
     call mat_create_block(C,Afull,3,2,3,5)

     sumElms = C%elms(23)+C%elms(24)+C%elms(25)+C%elms(28)+C%elms(29)+C%elms(30)

     call AssertEquals(sumElms,14.1d0,1.0d-12,"Mat_create_block not working properly")


!=================================
!    Unit test - mat_add_block:
!=================================

     Afull = 1.25d0
     !call mat_print(C, 1, C%nrow, 1, C%ncol, 6)
     call mat_add_block(C,Afull,3,2,2,6)
     !call mat_print(C, 1, C%nrow, 1, C%ncol, 6)

     sumElms = C%elms(27)+C%elms(28)+C%elms(29)+C%elms(32)+C%elms(33)+C%elms(34)

     call AssertEquals(sumElms,32.3d0,1.0d-12,"Mat_add_block not working properly")


!===================================
!    Unit test - mat_retrieve_block:
!===================================

     Afull = 0.0d0
     call mat_retrieve_block(C,Afull,3,2,2,6)

     sumElms = 0.0d0
     do i = 1, 3
        do j = 1, 2
           sumElms = sumElms + Afull(i,j)
        enddo
     enddo

     call AssertEquals(sumElms,32.3d0,1.0d-12,"Mat_retrieve_block not working properly")

     deallocate(Afull)

!===================================
!    Unit test - mat_zerohalf:
!===================================

     allocate(Ffull(7,7))

     call mat_zerohalf('lt',F)

     call mat_to_full(F,1.0d0,Ffull)

     sum1 = 0.0d0 ; sum2 = 0.0d0
     do i = 1, 7
        do j = 1, 7
           if (j > i) then
              sum1 = sum1 + Ffull(i,j)
           else
              sum2 = sum2 + Ffull(i,j)
           endif
        enddo
     enddo

     call AssertEquals(sum1,0.174866035580635d0,1.0d-12,"Mat_zerohalf/lower triangle not working properly")
     call AssertEquals(sum2,0.0d0,1.0d-12,"Mat_zerohalf/lower triangle not working properly")

     call mat_zerohalf('ut',F)
     call mat_zero(G)

     call AssertEquals(F,G,"Mat_zerohalf/upper triangle not working properly")

     deallocate(Ffull)

!====================================
!    Unit test - mat_report_sparsity:
!====================================

     !call mat_print(H, 1, H%nrow, 1, H%ncol, 6)
     call mat_report_sparsity(H,"test matrix",nnz,6)

     call AssertEquals(nnz,22,"Mat_report_sparsity not working properly")

!=====================================================
!    Unit test - mat_write_to_disk/mat_read_from_disk:
!=====================================================

     call mat_copy(1.0d0,H,F)

     open(1,file='testfile',status='new',action='write',form='unformatted',iostat=io)
     call AssertEquals(io,0,"Fortran OPEN not working properly")

     call mat_write_to_disk(1,H)

     close(1,status='keep',iostat=io)
     call AssertEquals(io,0,"Fortran CLOSE not working properly")

     open(2,file='testfile',status='old',action='read',form='unformatted',iostat=io)

     call mat_read_from_disk(2,G)

     close(2,status='delete')

     call AssertEquals(F,G,"Mat_write_to_disk/mat_read_from_disk not working properly")

!=====================================================
!    Unit test - mat_write_to_disk2/mat_read_from_disk2:
!=====================================================

     call mat_copy(1.0d0,H,F)

     open(1,file='testfile2',status='new',action='write')

     call mat_write_to_disk2(1,H)

     close(1,status='keep',iostat=io)

     open(2,file='testfile2',status='old',action='read')

     call mat_read_from_disk2(2,G)

     close(2,status='delete')

     call AssertEquals(F,G,"Mat_write_to_disk2/mat_read_from_disk2 not working properly")

!=================================
!    Unit test - mat_free:
!=================================

      call mat_free(A)
      !Test if elms pointer is disassociated:
      call AssertEquals(ASSOCIATED(A%elms),.false.,"Elms pointer not disassociated in mat_dense_free")

      call mat_free(B)
      call mat_free(C)
      call mat_free(D)
      call mat_free(E)
      call mat_free(F)
      call mat_free(G)
      call mat_free(H)
      call mat_free(L)
      call mat_free(M)
      call mat_free(N)
      call mat_free(O)
      call mat_free(P)
      call mat_free(Q)

!======================================
!    End of test suite - print status:
!======================================
      call print_status("dense matrix operations")

   end subroutine MatrixTestDense

!=================================================================
!    Helper routines for dense matrix unit tests
!=================================================================

!> \brief Simple matrix multiplication c = a*b
!> \author S. Host
!> \date 07-09-2010
subroutine my_mat_mul(a,b,c)
implicit none
     !> First argument
     type(matrix), intent(in) :: a
     !> Second argument
     type(matrix), intent(in) :: b
     !> Result - output
     type(matrix), intent(inout) :: c
     integer                  :: i, j, k
     real(realk)              :: sum

     call mat_zero(C)

     do i = 1, A%nrow
        do j = 1, B%ncol
           sum = 0
           do k = 1, A%ncol
              sum = sum + A%elms(A%nrow*(k-1)+i)*B%elms(B%nrow*(j-1)+k)
           enddo
           C%elms(C%nrow*(j-1)+i) = sum
        enddo
     enddo

end subroutine my_mat_mul

!> \brief Simple matrix addition c = a + b
!> \author S. Host
!> \date 07-09-2010
subroutine my_mat_add(a,b,c)
implicit none
     !> First argument
     type(matrix), intent(in) :: a
     !> Second argument
     type(matrix), intent(in) :: b
     !> Result - output
     type(matrix), intent(inout) :: c
     integer                  :: i

     call mat_zero(C)

     if (A%nrow /= B%nrow .or. A%nrow /= C%nrow) Stop 'Row dimensions wrong in my_mat_add'
     if (A%ncol /= B%ncol .or. A%ncol /= C%ncol) Stop 'Column dimensions wrong in my_mat_add'

     do i = 1, A%nrow*A%ncol
        C%elms(i) = A%elms(i) + B%elms(i)
     enddo

end subroutine my_mat_add

end module unittest_mat_dense

