!> @file 
!> Contains DIIS and EDIIS modules.

!> \brief DIIS (Direct Inversion in the Iterative Subspace) module. P. Pulay, CPL 73 (2), 393, 1980
!> \author L. Thogersen. Documented by S. Host.
!> \date 2003
module LINSCA_DIIS
  use av_utilities
  use scf_stats, only: stat_current_iteration, stat_tab
  use memory_handling
  use matrix_module
  use matrix_operations
  use precision
  private
  public :: DIIS, simple_averaging
  contains
!> \brief Construct average F and D from subspace from previous SCF iterations.
!> \author L. Thogersen
!> \date 2003
!> \param av Used to store info about SCF averaging
!> \param queue Subspace of Fock/KS and density matrices from previous SCF iterations
!> \param Dav The averaged density matrix
!> \param Fav The averaged Fock/KS matrix
  subroutine DIIS(av,queue,Dav,Fav)
    implicit none
    type(avItem), intent(inout)           :: av
    TYPE(util_HistoryStore),intent(in) :: queue
    type(Matrix) :: Dav,Fav
    real(realk),pointer  :: weights(:)
    real(realk) :: gdnorm,gdsqnorm, gdsnorm
    integer :: i, pos, msize

    msize = MIN(queue%used_entries,av%diis_history_size)
    pos = av%diis_pos
    if (queue%used_entries == 1) then
      call mat_assign(Dav,queue%D(1))
      call mat_assign(Fav,queue%F(1))
    else
!** Find DIIS weights
      call mem_alloc(weights,msize)
      call DIIS_WEIGHTS(av,msize,queue,weights)
      if (av%info_weight_final .or. av%info_diis) then
        WRITE(av%LUPRI,*) 'WEIGHTS,POS',POS
        call LS_OUTPUT(weights,1,msize,1,1,msize,1,1,av%lupri)
      endif
!** Construct average density and fock matrix
      if (.not. av%trilevel_gcscf_skip) then
         call get_AVERAGE_arr(av,'D',queue,av%diis_history_size,weights,Dav)
         call get_AVERAGE_arr(av,'F',queue,av%diis_history_size,weights,Fav)
      endif
      !call PURIFY(Ndim,Dav,S) !If new D is not found diagonalizing F
!** note in general info table that diis has been carried
!** out by putting a -1
      stat_tab(stat_current_iteration+1,6) = -1
      call mem_dealloc(weights)
    endif

  end subroutine DIIS

!> \brief Get weights for averaging Fock/KS and density matrices in DIIS.
!> \author L. Thogersen
!> \date 2003
!> \param av Used to store info about SCF averaging
!> \param msize Size of subspace
!> \param queue Subspace of Fock/KS and density matrices from previous SCF iterations
!> \param weights The weights for averaging Fock/KS and density matrices
  subroutine DIIS_WEIGHTS(av,msize,queue,weights)
    implicit none
    type(avItem), intent(inout) :: av
    integer, intent(in) :: msize
    type(util_historyStore)  :: queue
    real(realk), intent(out) :: weights(msize)
    real(realk) :: Amat(msize+1,msize+1)
    real(realk) :: bvec(msize+1)
    integer :: i,ierr,ipiv(msize+1)
    ierr=0

    Amat(1:msize,1:msize) = queue%DIIS_Amat(1:msize,1:msize)
    do i = 1, msize
      Amat(msize+1,i) = -1E0_realk
      Amat(i,msize+1) = -1E0_realk
    end do
    Amat(msize+1,msize+1) = 0.0
    if (av%info_diis) then
      write(av%lupri,*)'DIIS A matrix'
      call ls_output(Amat,1,msize+1,1,msize+1,msize+1,msize+1,1,av%lupri)
    endif    
    bvec(1:msize) = 0.0E0_realk
    bvec(msize+1) = -1.0E0_realk
    call DGESV(msize+1,1,Amat,msize+1,ipiv,bvec,msize+1,ierr)
    if (ierr == 0) then
      weights(1:msize) = bvec(1:msize)
    else
      if (av%trilevel_gcscf) then
        av%trilevel_gcscf_skip = .true.
        write(av%lupri,*) 'GCSCF mode: DIIS failed and is skipped'
      else
        WRITE(av%LUPRI,*) 'Error ',ierr,'  in diis DGESV'
        STOP 'error in DIIS DGESV'
      endif
    endif

  end subroutine DIIS_WEIGHTS

!> \brief Simple averaging of F and D: Fav = 1/2*(F(n) + F(n-1)), used for Van Lenthe Scheme
!> \author S. Host
!> \date march 2010
!> \param av Used to store info about SCF averaging
!> \param iteration Current SCF iteration
!> \param Dav The averaged density matrix
!> \param Fav The averaged Fock/KS matrix
  subroutine simple_averaging(av,iteration,Dav,Fav)
    implicit none
    type(avItem), intent(in)    :: av
    integer,intent(in)          :: iteration
    type(Matrix), intent(inout) :: Dav,Fav

    if (iteration == 1) then
       !nothing
    else
       ! Fnew = 0.5 (Fnew + Fold), Dnew = 0.5 (Dnew + Dold)
       call mat_scal(0.5E0_realk, Fav)
       call mat_daxpy(0.5E0_realk, av%Fprev, Fav) 
       call mat_scal(0.5E0_realk, Dav)
       call mat_daxpy(0.5E0_realk, av%Dprev, Dav)
    endif

  end subroutine simple_averaging

end module LINSCA_DIIS

! commentet out by TK - do not seem to work and no testcase. TK
!!$!> EDIIS (energy-DIIS) module. K. N. Kudin, G. E. Scuseria et al., JCP 116, 8255 (2002)
!!$!> \author L. Thogersen
!!$!> \date 2003
!!$module LINSCA_EDIIS
!!$  use av_utilities
!!$  use memory_handling
!!$  use matrix_module
!!$  use matrix_operations
!!$  use precision
!!$  contains
!!$!> \brief Modified DIIS. Averaging is based on energy minimization.
!!$!> \author L. Thogersen
!!$!> \date 2003
!!$!> \param av Used to store info about SCF averaging
!!$!> \param queue Subspace of Fock/KS and density matrices from previous SCF iterations
!!$!> \param Dav The averaged density matrix
!!$!> \param Fav The averaged Fock/KS matrix
!!$  subroutine EDIIS(av,queue,Dav,Fav)
!!$    use scf_stats, only: stat_current_iteration, stat_tab
!!$    implicit none
!!$    type(avItem), intent(inout)           :: av
!!$    TYPE(util_HistoryStore),intent(in) :: queue
!!$    type(Matrix) :: Dav,Fav
!!$    integer :: msize
!!$    real(realk), pointer :: weights(:)
!!$
!!$    msize = MIN(queue%used_entries,av%ediis_history_size)
!!$    if (queue%used_entries == 1) then
!!$      call mat_assign(Dav,queue%D(1))
!!$      call mat_assign(Fav,queue%F(1))
!!$    else
!!$!** Find EDIIS weights
!!$      call mem_alloc(weights,msize)
!!$      call EDIIS_WEIGHTS(av,msize,queue,weights)
!!$      if (av%info_weight_final .or. av%info_ediis) then
!!$        WRITE(av%LUPRI,*) 'WEIGHTS,POS',av%ediis_POS
!!$        call LS_OUTPUT(weights,1,msize,1,1,msize,1,1,av%lupri)
!!$      endif
!!$!** Construct average density and fock matrix
!!$      call get_AVERAGE_arr(av,'D',queue,av%ediis_history_size,weights,Dav)
!!$      call get_AVERAGE_arr(av,'F',queue,av%ediis_history_size,weights,Fav)
!!$      !call PURIFY(Ndim,Dav,S) !If new D is not found diagonalizing F
!!$!** note in general info table that ediis has been carried
!!$!** out by putting a -2
!!$      stat_tab(stat_current_iteration+1,6) = -2
!!$      call mem_dealloc(weights)
!!$    endif
!!$
!!$  end subroutine EDIIS
!!$
!!$!> \brief Get weights for averaging Fock/KS and density matrices in EDIIS.
!!$!> \author L. Thogersen
!!$!> \date 2003
!!$!> \param av Used to store info about SCF averaging
!!$!> \param msize Size of subspace
!!$!> \param queue Subspace of Fock/KS and density matrices from previous SCF iterations
!!$!> \param weights The weights for averaging Fock/KS and density matrices
!!$  subroutine EDIIS_WEIGHTS(av,msize,queue,weights)
!!$    implicit none
!!$    type(avItem), intent(in)           :: av
!!$    integer, intent(in) :: msize
!!$    TYPE(util_HistoryStore),intent(in) :: queue
!!$    real(realk), intent(out) :: weights(msize)
!!$    real(realk) :: grad(msize-1), hes(msize-1,msize-1), Amat(msize-1,msize-1), bvec(msize-1)
!!$    type(Matrix) :: Dav, Fav
!!$    integer :: ipiv(msize-1), error, i, j, his_start
!!$
!!$    his_start = queue%current_position
!!$!** Find Amat and bvec
!!$    call get_A_and_b(av,msize,queue,Amat,bvec)
!!$!** Find delta c
!!$    !** Solve system of linear equations - solution vector put in bvec
!!$    CALL DGESV(msize-1,1,Amat,msize-1,ipiv,bvec,msize-1,error)  !LAPACK
!!$    if (error /= 0) then 
!!$       write(av%lupri,*) &
!!$            &'  EDIIS_WEIGHTS: Error', error
!!$       call ls_output(Amat,1,msize-1,1,msize-1,msize-1,msize-1,1,av%lupri)
!!$       stop 'EDIIS_WEIGHTS: error in EDIIS'
!!$    end if
!!$    if (av%info_ediis) then
!!$       write(av%lupri,*)'Solution vector'
!!$       call ls_output(bvec,1,msize-1,1,1,msize-1,1,1,av%lupri)
!!$    endif
!!$    j = 0
!!$    do i = 1,msize
!!$      if (i == av%ediis_pos) then
!!$        weights(i) = 1.0E0_realk - SUM(bvec)
!!$      else
!!$        j = j + 1
!!$        weights(i) = bvec(j)
!!$      endif
!!$    enddo   
!!$
!!$  end subroutine EDIIS_WEIGHTS

!!$!> \brief Construct A matrix and b vector for EDIIS (see paper).
!!$!> \author L. Thogersen
!!$!> \date 2003
!!$!> \param av Used to store info about SCF averaging
!!$!> \param msize Size of subspace
!!$!> \param queue Subspace of Fock/KS and density matrices from previous SCF iterations
!!$!> \param Amat The A matrix
!!$!> \param bvec The b vector
!!$  subroutine get_A_and_b(av,msize,queue,Amat,bvec)
!!$    implicit none
!!$    type(avItem), intent(in) :: av
!!$    integer, intent(in)      :: msize
!!$    type(util_HistoryStore)  :: queue
!!$    real(realk), intent(out) :: Amat(msize-1,msize-1),bvec(msize-1)
!!$    integer :: his_start,i,j,k,m,i2,j2
!!$    real(realk) :: TrDoFo, TrDjFo,TrDoFj,coef(msize),grad(msize),hes(msize,msize),TrDiFj
!!$
!!$    his_start = queue%current_position
!!$    coef = 0.0E0_realk
!!$    coef(av%ediis_pos) = 1E0_realk
!!$    !if (av%debug_ediis) then
!!$    !  call debug_EDbar(msize,coef,his_start,queue,grad,hes)
!!$    !  k = 0
!!$    !  do i = 1,msize
!!$    !    if (i /= av%ediis_pos) then
!!$    !      k = k+1
!!$    !      bvec(k) = grad(i) - grad(av%ediis_pos) 
!!$    !      m = 0
!!$    !      do j = 1,msize
!!$    !        if (j /= av%ediis_pos) then
!!$    !          m = m+1
!!$    !          Amat(k,m) = hes(i,j) - hes(i,av%ediis_pos) - hes(av%ediis_pos,j) + hes(av%ediis_pos,av%ediis_pos)
!!$    !        endif
!!$    !      enddo
!!$    !    endif
!!$    !  enddo
!!$    !  WRITE(av%LUPRI,*) 'Finite difference derivatives'
!!$    !  write(av%lupri,*) 'EDIIS gradient '
!!$    !  call ls_output(bvec,1,msize-1,1,1,msize-1,1,1,av%lupri)
!!$    !  write(av%lupri,*) 'EDIIS hessian '
!!$    !  call ls_output(Amat,1,msize-1,1,msize-1,msize-1,msize-1,1,av%lupri)
!!$    !endif
!!$
!!$!** Construct bvec
!!$    !   b_i = 2TrDoFo - 2TrD_i Fo
!!$    j = queue%current_position
!!$    i = av%ediis_pos
!!$    TrDoFo = mat_dotproduct(queue%D(j),queue%F(j))
!!$    do k = 1,msize-1
!!$      j = j-1; if(j==0) j=av%cfg_settings(av%cfg_set_type)%max_history_size
!!$      i = i-1; if(i==0) i=av%ediis_history_size-1
!!$      bvec(i) = 2.0E0_realk*(TrDoFo - mat_dotproduct(queue%D(j),queue%F(his_start)))
!!$    enddo
!!$!** Construct Amat
!!$    !   A_ij = Tr(Di-Do)(Fj-Fo) + Tr(Dj-Do)(Fi-Fo) = 2TrDoFo + TrDiFj + TrDjFi - TrDiFo - TrDjFo - TrDoFj - TrDoFi
!!$    j = queue%current_position
!!$    i = av%ediis_pos
!!$    Amat = 2E0_realk*TrDoFo  ! 2TrDoFo
!!$    do k = 1,msize-1  
!!$      j = j-1; if(j==0) j=av%cfg_settings(av%cfg_set_type)%max_history_size
!!$      i = i-1; if(i==0) i=av%ediis_history_size-1
!!$      TrDoFj = mat_dotproduct(queue%D(his_start),queue%F(j)) 
!!$      Amat(:,i) = Amat(:,i) - TrDoFj                                         ! -TrDoFj
!!$      Amat(i,:) = Amat(i,:) - TrDoFj                                         ! -TrDoFi
!!$      TrDjFo = mat_dotproduct(queue%D(j),queue%F(his_start))
!!$      Amat(i,:) = Amat(i,:) - TrDjFo                                         ! -TrDiFo
!!$      Amat(:,i) = Amat(:,i) - TrDjFo                                         ! -TrDjFo
!!$      j2 = queue%current_position
!!$      i2 = av%ediis_pos
!!$      do m = 1,msize-1
!!$        j2 = j2-1; if(j2==0) j2=av%cfg_settings(av%cfg_set_type)%max_history_size
!!$        i2 = i2-1; if(i2==0) i2=av%ediis_history_size-1
!!$        TrDiFj = mat_dotproduct(queue%D(j2),queue%F(j))
!!$        Amat(i2,i) = Amat(i2,i) + TrDiFj                                     ! TrDiFj
!!$        Amat(i,i2) = Amat(i,i2) + TrDiFj                                     ! TrDjFi
!!$      enddo
!!$    enddo
!!$    if (av%debug_ediis) then
!!$      WRITE(av%LUPRI,*) 'EDIIS derivatives'
!!$      write(av%lupri,*) 'gradient '
!!$      call ls_output(bvec,1,msize-1,1,1,msize-1,1,1,av%lupri)
!!$      write(av%lupri,*) 'hessian '
!!$      call ls_output(Amat,1,msize-1,1,msize-1,msize-1,msize-1,1,av%lupri)
!!$    endif       
!!$
!!$  end subroutine get_A_and_b
!!$
!!$end module LINSCA_EDIIS

