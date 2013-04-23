MODULE Fock_evaluator
  use precision
  use matrix_module
  use matrix_operations
  use dal_interface, only: di_get_fock_LSDALTON
  use lsdalton_fock_module
  use TYPEDEFTYPE, only: lssetting, lsitem

! This module evaluates the Fock/Kohn-Sham matrix using chosen
! algorithm, matrix representation, etc.
!
!   PUBLIC

CONTAINS

   SUBROUTINE FCK_get_fock(d, F, Etotal)
      IMPLICIT NONE
      TYPE(Matrix), intent(in)    :: D
      type(matrix), intent(inout) :: F  !output 
      real(realk), INTENT(OUT) :: Etotal
      type(matrix) :: Dtmp(1),Fock(1)
      real(realk) :: E(1)
      integer :: ndmat
      ndmat = 1
      call mat_init(Dtmp(1),D%nrow,D%ncol)
      call mat_init(Fock(1),F%nrow,F%ncol)
      call mat_assign(Dtmp(1),D)
      call mat_assign(Fock(1),F)
      !cfg_nfock = cfg_nfock + 1
      CALL di_get_fock_LSDALTON(Dtmp,lsint_fock_data%H1,Fock,ndmat,E, &
           & lsint_fock_data%lupri,lsint_fock_data%luerr,  &
           & lsint_fock_data%ls)
      Etotal = E(1)
      call mat_assign(F,Fock(1))
      call mat_free(Dtmp(1))
      call mat_free(Fock(1))
   END SUBROUTINE FCK_get_fock

   subroutine fck_scale_virt_fock(S,H1,F,D)
     !to create a N-1 potential for the virtual orbitals
     !allowing for more 'occupied-like' virtual orbitals
     !good when configuration shifts are needed!!!
     !** F_mod = F + (n-1/n - 1)Q^T (F - H1) Q
     type(Matrix), intent(in) :: S,H1,D
     type(Matrix), intent(inout) :: F
     type(matrix) :: Q,wrk1,wrk2
     real(realk) :: konst
     integer :: ndim
     call lsquit('fck_scale_virt_fock not tested - so now I am sure it is wrong',-1)
     !this subroutine is not tested - so now I am sure it is wrong
     !it uses cfg_nocc before it is set - if it is ever set. TK 2012
!!$     ndim = S%nrow
!!$     call mat_init(Q,ndim,ndim)
!!$     call mat_init(wrk1,ndim,ndim)
!!$     call mat_init(wrk2,ndim,ndim)
!!$
!!$     !** wrk2 = F - h1
!!$     call mat_add(1E0_realk,F,-1E0_realk,H1,wrk2)
!!$     !** Q = 1 - DS
!!$     call mat_mul(D,S,'n','n',1E0_realk,0E0_realk,wrk1)
!!$     call mat_add_identity(1E0_realk,-1E0_realk,wrk1,Q)
!!$     !** wrk1 = (F - h1) Q = wrk2 Q
!!$     call mat_mul(wrk2,Q,'n','n',1E0_realk,0E0_realk,wrk1)
!!$     !** konst = -1/n
!!$     konst = -1E0_realk/(2E0_realk*cfg_nocc)
!!$     !** F_mod = F + konst Q^T (F - H1) Q = F + konst Q^T wrk1
!!$     call mat_mul(Q,wrk1,'T','n',konst,1E0_realk,F)
!!$
!!$     call mat_free(Q)
!!$     call mat_free(wrk1)
!!$     call mat_free(wrk2)
   end subroutine fck_scale_virt_fock

   ! If the virtual orbitals are scaled in the calculation
   ! then this routine unscales them such that the "actual" 
   ! Fock matrix is outputted.
   subroutine fck_unscale_virt(H1,S,D,F,nocc)
     ! F = F_mod - k Q^T G(D) Q = F_mod - k/(1+k)Q^T (F_mod - h1) Q  where
     ! k = (n-1)/n - 1 = -1/n
     implicit none
     type(matrix), intent(in)    :: H1,S,D
     type(Matrix),intent(inout)  :: F
     !> Number of occupied orbitals
     integer, intent(in) :: nocc
     type(Matrix) :: Q,wrk1,wrk2
     real(realk) :: konst,Etotal
     integer :: ndim

     if (.false.) then
     ndim = S%nrow
     call mat_init(Q,ndim,ndim)
     call mat_init(wrk1,ndim,ndim)
     call mat_init(wrk2,ndim,ndim)
     !wrk2 = F_mod - h1
     call mat_add(1E0_realk,F,-1E0_realk,H1,wrk2)
     !Q = 1 - DS
     call mat_mul(D,S,'n','n',1E0_realk,0E0_realk,wrk1)
     call mat_add_identity(1E0_realk,-1.0E0_realk,wrk1,Q)
     !wrk1 = (F_mod - h1) Q = wrk2 Q
     call mat_mul(wrk2,Q,'n','n',1E0_realk,0E0_realk,wrk1)
     !F = F_mod + konst Q^T (F_mod - h1) Q = F_mod + konst Q^T wrk1
     konst = 1E0_realk/(nocc*2E0_realk - 1E0_realk)
     call mat_mul(Q,wrk1,'t','n',konst,1E0_realk,F)
     call mat_free(Q)
     call mat_free(wrk1)
     call mat_free(wrk2)
     else
       !cfg_scale_virt = .false.
       CALL lsquit('replaced di_get_fock(d,F,Etotal) with this quit statement',-1)
     endif
   END SUBROUTINE fck_unscale_virt
  
END MODULE Fock_evaluator
