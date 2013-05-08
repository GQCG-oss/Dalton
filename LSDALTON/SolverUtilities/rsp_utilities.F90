MODULE rsp_util
! This module contains stuff that has to be separate because F90 allows
! only tree-type dependencies between files and there is no way to resolve
! circular dependency by eg. creation of a header file.
   USE Matrix_operations
   USE Matrix_operations_aux ! mat_mo_precond,mat_new_mo_precond,mat_new_complex_precond,..
   private :: Cmo_final,  Cmo_final_saved, orbE_final
   public
   LOGICAL, SAVE :: flushed, Cmo_final_saved
   type(Matrix), save :: Cmo_final
   real(realk), pointer, save :: orbE_final(:)

CONTAINS
  subroutine init_rsp_util()
    implicit none 
    flushed = .false.
    Cmo_final_saved = .false.
  end subroutine init_rsp_util

!************************************
!** PRECONDITIONING
!************************************
    subroutine util_save_MOinfo(F,S,nocc)
      implicit none 
      type(Matrix), intent(in) :: F,S
      integer :: i,j,ndim,nocc
      Cmo_final_saved = .true.
      ndim=F%nrow
      call mat_init(Cmo_final,ndim,ndim)
      ALLOCATE(orbE_final(S%nrow*2))!allow for unrestricted
      call mat_diag_f(F,S,orbE_final, Cmo_final)

      !fixme: uhf case needs more work?

      WRITE(6,*) 
      WRITE(6,'("*******************************************************************")')
      WRITE(6,'(" Orbital energies     ",i3," occupied orbitals")') nocc
      j = 1
      do
         if (j+4 > nocc) exit 
         WRITE(6,"(5f13.7)") (orbE_final(i),i=j,j+4)
         j = j+5
      enddo
      WRITE(6,"(5f13.7)") (orbE_final(i),i=j,nocc)
      
      WRITE(6,*)
      WRITE(6,'("                      ",i3," virtual orbitals")') ndim - nocc
      j = nocc + 1
      do 
         if (j+4 > ndim) exit
         WRITE(6,"(5f13.7)") (orbE_final(i),i=j,j+4)
         j = j+5
      enddo
      WRITE(6,"(5f13.7)") (orbE_final(i),i=j,ndim)
      WRITE(6,*)
      ! we keep Cmo_final and orbE_final in the memory!
    end subroutine util_save_MOinfo

    subroutine util_get_CMO(CMO)
      implicit none 
      type(Matrix), intent(inout) :: CMO
      IF(.NOT.Cmo_final_saved)call lsquit('Cmo_final have not been saved (util_save_MOinfo)',-1)
      call mat_assign(CMO,Cmo_final)
    end subroutine util_get_CMO

    subroutine get_rsp_trials_from_MO(nroots,bvec_ao,nocc)
      implicit none
      integer, intent(in) :: nroots,nocc
      type(Matrix), intent(inout) :: bvec_ao(nroots)
      type(matrix) :: tmp_mat
      real(realk) :: omegas(NOCC*(Cmo_final%nrow-NOCC)),x
      integer :: i,j,a,n,Mindex((Cmo_final%nrow-NOCC)*NOCC,2),xv(2),ndim
      real(realk), allocatable   :: iniguess_full(:,:) !FIXME: mat_create_elm does not work with BSM
      IF(.NOT.Cmo_final_saved)call lsquit('Cmo_final have not been saved (util_save_MOinfo)',-1)
      ndim = Cmo_final%nrow
      allocate(iniguess_full(ndim,ndim))
      IF(NROOTS.GT.(Cmo_final%nrow-NOCC)*NOCC)THEN
       Call lsquit("Error in get_rsp_trials_from_MO: (Cmo_final%nrow-NOCC)*NOCC .LT. NEXCIT",-1)
      ENDIF
      write(6,*) 'get_rsp_trials_from_MO, ndim = ', ndim,'nocc = ',nocc
      n = 0
      !Find the DeltaE_ai orbital energy differences as a first guess of 
      !the size of the excitation-energies 
      do i = 1,NOCC
        do a = NOCC+1,ndim
          n = n+1
          omegas(n) = orbE_final(a) - orbE_final(i)
          Mindex(n,1) = a
          Mindex(n,2) = i
        enddo
      enddo
      !sort them such that the lowest is first
      !SONIA FIXME: use better sorting algorithm
      do i = 1,n
        do j = 1,n-1
          if (omegas(j) > omegas(j+1)) then
            x = omegas(j+1)
            omegas(j+1) = omegas(j)
            omegas(j) = x
            xv = Mindex(j+1,:)
            Mindex(j+1,:) = Mindex(j,:)
            Mindex(j,:) = xv
          endif
        enddo
      enddo
      !SONIA: to be removed
      !if (info_rsp) write(6,*) 'Inside 1st trial: resorted deltaepsi'
      !do i=1,n
      !   write(6,*) 'I=',I, ' omega_I', omegas(i)
      !end do
      !SONIA: end of to be removed

      !choose the nroot lowest for trials
      !and create trial vectors
      call mat_init(tmp_mat,ndim,ndim)
      do i = 1,nroots
        !call mat_zero(Bvec_ao(i))
        !call mat_create_elm(Mindex(i,1),Mindex(i,2),1.0E0_realk,Bvec_ao(i))
        iniguess_full = 0.0E0_realk
        iniguess_full(Mindex(i,1),Mindex(i,2)) = 1.0E0_realk
        call mat_set_from_full(iniguess_full,1.0E0_realk, Bvec_ao(i), 'MO iniguess')
        !transform to AO-basis
        call mat_mul(Cmo_final,Bvec_ao(i),'n','n',1E0_realk,0E0_realk,tmp_mat)
        call mat_mul(tmp_mat,Cmo_final,'n','t',1E0_realk,0E0_realk,Bvec_ao(i))
      enddo
      call mat_free(tmp_mat)
      deallocate(iniguess_full)

    end subroutine get_rsp_trials_from_MO

    subroutine MO_precond(Gn,S,omega,Gnt,nocc)
      !Preconditioning of Gn by solving A Gnt = Gn
      !First transforming Gn to the MO-basis
      !Second divide all ai elements by 2(E_a - E_i)
      !Third transform the solution Gnt back to the AO-basis
      implicit none
      type(Matrix), intent(in) :: Gn,S
      !integer, intent(in) :: symm
      real(realk), intent(in) :: omega
      type(Matrix), intent(inout) :: Gnt !OUTPUT
      type(Matrix) :: Gn_MO
      integer :: ndim,nocc,nvirt
      logical :: cov

      ndim = Gn%nrow
      nvirt = ndim - nocc
!**** Transform Gn to the MO-basis. 
      call mat_init(Gn_MO,ndim,ndim)
      !C^T Gn C
      cov = .true.
      call util_AO_to_MO(S,Gn,Gn_MO,cov)
      !We are only interested in the ia part of the matrix
!**** Divide all ai elements by E_a - E_i
      call mat_mo_precond(nocc,omega,orbE_final,Gn_MO)
!**** Transform Gn_MO back to the AO-basis
      ! C Gn_MO C^T
      cov = .false.
      call util_MO_to_AO(S,Gn_MO,Gnt,cov)
      !Free matrices
      call mat_free(Gn_MO)
    end subroutine MO_precond

    subroutine new_std_MO_precond(rp,rm,S,omega,xp,xm,nocc)
      !Preconditioning of Gn by solving A Gnt = Gn
      !First transforming Gn to the MO-basis
      !Second divide all ai elements by 2(E_a - E_i)
      !Third transform the solution Gnt back to the AO-basis
      implicit none
      type(Matrix), intent(in) :: rp,rm,S
      !integer, intent(in) :: symm
      real(realk), intent(in) :: omega
      type(Matrix), intent(inout) :: xp,xm !OUTPUT
      type(Matrix) :: rp_MO,rm_mo
      integer :: ndim,nocc,nvirt
      logical :: cov

      ndim = xp%nrow
      nvirt = ndim - nocc
!**** Transform Gn to the MO-basis. 
      call mat_init(rp_MO,ndim,ndim)
      call mat_init(rm_MO,ndim,ndim)
      !C^T Gn C
      cov = .true.
      call util_AO_to_MO(S,rp,rp_MO,cov)
      call util_AO_to_MO(S,rm,rm_MO,cov)
      !We are only interested in the ia part of the matrix
!**** Divide all ai elements by E_a - E_i
      call mat_new_mo_precond(nocc,omega,orbE_final,rp_MO,rm_MO)
!**** Transform Gn_MO back to the AO-basis
      ! C Gn_MO C^T
      cov = .false.
      call util_MO_to_AO(S,rp_MO,xp,cov)
      call util_MO_to_AO(S,rm_MO,xm,cov)
      !Free matrices
      call mat_free(rp_MO)
      call mat_free(rm_MO)
    end subroutine new_std_MO_precond
    
    subroutine new_complex_precond(rp,rm,rpi,rmi,S,omega,gammma,xp,xm,xpi,xmi,nocc)
      !Preconditioning of Gn by solving A Gnt = Gn
      !First transforming Gn to the MO-basis
      !Second divide all ai elements by 2(E_a - E_i)
      !Third transform the solution Gnt back to the AO-basis
      implicit none
      type(Matrix), intent(in) :: rp,rm,rpi,rmi,S
      !integer, intent(in) :: symm
      real(realk), intent(in) :: omega,gammma
      type(Matrix), intent(inout) :: xp,xm,xpi,xmi !OUTPUT
      type(Matrix) :: rp_MO,rm_mo,rpi_mo,rmi_mo
      integer :: ndim,nocc,nvirt
      logical :: cov

      ndim = xp%nrow
      nvirt = ndim - nocc
!**** Transform Gn to the MO-basis. 
      call mat_init(rp_MO,ndim,ndim)
      call mat_init(rm_MO,ndim,ndim)
      call mat_init(rpi_MO,ndim,ndim)
      call mat_init(rmi_MO,ndim,ndim)
      !C^T Gn C
      cov = .true.
      call util_AO_to_MO(S,rp,rp_MO,cov)
      call util_AO_to_MO(S,rm,rm_MO,cov)
      call util_AO_to_MO(S,rpi,rpi_MO,cov)
      call util_AO_to_MO(S,rmi,rmi_MO,cov)
      !We are only interested in the ia part of the matrix
!**** Divide all ai elements by E_a - E_i
      call mat_new_complex_precond(nocc,omega,gammma,orbE_final,rp_MO,rm_MO,rpi_mo,rmi_mo)
!**** Transform Gn_MO back to the AO-basis
      ! C Gn_MO C^T
      cov = .false.
      call util_MO_to_AO(S,rp_MO,xp,cov)
      call util_MO_to_AO(S,rm_MO,xm,cov)
      call util_MO_to_AO(S,rpi_MO,xpi,cov)
      call util_MO_to_AO(S,rmi_MO,xmi,cov)
      !Free matrices
      call mat_free(rp_MO)
      call mat_free(rm_MO)
      call mat_free(rpi_MO)
      call mat_free(rmi_MO)
    end subroutine new_complex_precond
   
 subroutine MO_precond_complex(Gn,D,S,Gnt,nocc)
      !Preconditioning of Gn by solving A Gnt = Gn
      !First transforming Gn to the MO-basis
      !Second divide all ai elements by 2(E_a - E_i)
      !Third transform the solution Gnt back to the AO-basis
      implicit none
      type(Matrix), intent(in) :: Gn,S,D
      !integer, intent(in) :: symm
      type(Matrix), intent(inout) :: Gnt !OUTPUT
      type(Matrix) :: Gn_MO
      integer :: ndim,nocc,nvirt
      logical :: cov

      ndim = Gn%nrow
      nvirt = ndim - nocc
!**** Transform Gn to the MO-basis. 
      call mat_init(Gn_MO,ndim,ndim)
      !C^T Gn C
      cov = .true.
      call util_AO_to_MO(S,Gn,Gn_MO,cov)
      !We are only interested in the ia part of the matrix
!**** Divide all ai elements by E_a - E_i
      call mat_mo_precond_complex(nocc,orbE_final,Gn_MO)
!**** Transform Gn_MO back to the AO-basis
      ! C Gn_MO C^T
      cov = .false.
      call util_MO_to_AO(S,Gn_MO,Gnt,cov)
      !Free matrices
      call mat_free(Gn_MO)
      call util_scriptPx('N',D,S,Gnt)
      end subroutine MO_precond_complex

    subroutine util_AO_to_MO(S,X_ao,X_mo,cov)
      ! if covariant
      !Transform from the AO covariant basis to the MO basis
      !X_mo = C^T X_ao C
      ! else contravariant
      !Transform from the AO contravariant basis to the MO basis
      !X_mo = C^T S X_ao S C
      implicit none
      logical :: cov
      type(Matrix), intent(in) :: S,X_ao
      type(Matrix), intent(inout) :: X_mo  !output
      type(Matrix) :: X,SC
      integer :: ndim
      IF(.NOT.Cmo_final_saved)call lsquit('Cmo_final have not been saved (util_save_MOinfo)',-1)

      call mat_init(X,S%nrow,S%ncol)
      call mat_init(SC,S%nrow,S%ncol)
      ndim = S%nrow
      if (cov) then
        call mat_assign(SC,Cmo_final)
      else
        call mat_mul(S,Cmo_final,'n','n',1E0_realk,0E0_realk,SC)
      endif
      call mat_mul(X_ao,SC,'n','n',1E0_realk,0E0_realk,X)
      call mat_mul(SC,X,'t','n',1E0_realk,0E0_realk,X_mo)

      call mat_free(X)
      call mat_free(SC)

    end subroutine util_AO_to_MO

    subroutine util_MO_to_AO(S,X_mo,X_ao,cov)
      ! if covariant
      !Transform from MO to covariant AO basis
      !  X_ao =  S C X_mo C^T S
      ! else contravariant
      !Transform from MO to contravariant AO basis
      !X_ao = C X_mo C^T
      implicit none
      logical :: cov
      type(Matrix), intent(in) :: S,X_mo
      type(Matrix), intent(inout) :: X_ao  !output
      type(Matrix) :: X, SC
      IF(.NOT.Cmo_final_saved)call lsquit('Cmo_final have not been saved (util_save_MOinfo)',-1)
      call mat_init(X,X_mo%nrow,X_mo%ncol)
      call mat_init(SC,S%nrow,S%ncol)

      if (cov) then
        call mat_mul(S,Cmo_final,'n','n',1E0_realk,0E0_realk,SC)
      else
        call mat_assign(SC,Cmo_final)
      endif
      call mat_mul(X_MO,SC,'n','t',1E0_realk,0E0_realk,X)
      call mat_mul(SC,X,'n','n',1E0_realk,0E0_realk,X_ao)

      call mat_free(X)
      call mat_free(SC)
      
    end subroutine util_MO_to_AO

    subroutine util_free_mostuff()
      implicit none
      IF(.NOT.Cmo_final_saved)call lsquit('Cmo_final have not been saved (util_save_MOinfo)',-1)
      call mat_free(cmo_final)
      deallocate(orbE_final)
      Cmo_final_saved = .FALSE.
    end subroutine util_free_mostuff
!*******************
!** Projecting PXQ^T + QXP^T
!*******************************
!If trans = T or t scriptPT is made:
!scriptP^T(X) = P^T X Q + Q^T X P in orbitalbase
!If trans = N or n scriptP is made:
!scriptP(X) = P X Q^T + Q X P^T in orbitalbase
    subroutine util_scriptPx(trans,D,S,X)
      implicit none
      character, intent(in) :: trans
      type(Matrix), intent(in) :: D,S
      type(Matrix), intent(inout) :: x
      double precision :: norm
      type(Matrix) :: P, Q, scr1,scr2,scr3
    
      call mat_init(P,d%nrow,D%ncol)
      call mat_init(Q,D%nrow,D%ncol)
      call mat_init(scr1,D%nrow,D%ncol)
      call mat_init(scr2,D%nrow,D%ncol)
      call mat_init(scr3,D%nrow,D%ncol)
      call mat_mul(D,S,'n','n',1E0_realk,0E0_realk,P)    
      call mat_add_identity(1.0E0_realk,-1.0E0_realk,P,Q)  !Q = I-P
      if (trans == 'T' .or. trans == 't') then
        !** scr2 = P^TXQ
        call mat_mul(P,x,'t','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,Q,'n','n',1E0_realk,0E0_realk,scr2)
        !** scr3 = Q^TXP
        call mat_mul(Q,X,'t','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,P,'n','n',1E0_realk,0E0_realk,scr3)
        !** X = P^TXQ + Q^TXP
        call mat_add(1E0_realk,scr2,1E0_realk,scr3,X)
      elseif (trans == 'N' .or. trans == 'n') then
        !** scr2 = PXQ^T
        call mat_mul(P,x,'n','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,Q,'n','t',1E0_realk,0E0_realk,scr2)
        !** scr3 = QXP^T
        call mat_mul(Q,X,'n','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,P,'n','t',1E0_realk,0E0_realk,scr3)
        !** X = PXQ^T + QXP^T
        call mat_add(1E0_realk,scr2,1E0_realk,scr3,X)
      else
        WRITE(6,*) 'unknown transformation in scriptPx',trans
        STOP 'unknown transformation in scriptPx'
      endif
      call mat_free(P)
      call mat_free(Q)
      call mat_free(scr1)
      call mat_free(scr2)
      call mat_free(scr3)
    
    end subroutine util_scriptPx

    subroutine util_scriptPxOAO(trans,D,X)
      implicit none
      character, intent(in)       :: trans
      type(Matrix), intent(in)    :: D
      type(Matrix), intent(inout) :: x
      double precision            :: norm
      type(Matrix)                :: Q, scr1,scr2,scr3
    
      call mat_init(Q,D%nrow,D%ncol)
      call mat_init(scr1,D%nrow,D%ncol)
      call mat_init(scr2,D%nrow,D%ncol)
      call mat_init(scr3,D%nrow,D%ncol)
      call mat_add_identity(1.0E0_realk,-1.0E0_realk,D,Q)  !Q = I-P
      if (trans == 'T' .or. trans == 't') then
        !** scr2 = P^TXQ
        call mat_mul(D,x,'t','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,Q,'n','n',1E0_realk,0E0_realk,scr2)
        !** scr3 = Q^TXP
        call mat_mul(Q,X,'t','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,D,'n','n',1E0_realk,0E0_realk,scr3)
        !** X = P^TXQ + Q^TXP
        call mat_add(1E0_realk,scr2,1E0_realk,scr3,X)
      elseif (trans == 'N' .or. trans == 'n') then
        !** scr2 = PXQ^T
        call mat_mul(D,x,'n','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,Q,'n','t',1E0_realk,0E0_realk,scr2)
        !** scr3 = QXP^T
        call mat_mul(Q,X,'n','n',1E0_realk,0E0_realk,scr1)
        call mat_mul(scr1,D,'n','t',1E0_realk,0E0_realk,scr3)
        !** X = PXQ^T + QXP^T
        call mat_add(1E0_realk,scr2,1E0_realk,scr3,X)
      else
        WRITE(6,*) 'unknown transformation in scriptPx',trans
        STOP 'unknown transformation in scriptPx'
      endif
      call mat_free(Q)
      call mat_free(scr1)
      call mat_free(scr2)
      call mat_free(scr3)
    
    end subroutine util_scriptPxOAO
        
END MODULE rsp_util

